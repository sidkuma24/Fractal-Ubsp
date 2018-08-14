#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct list_node {
    unsigned short qalfa;
    unsigned short qbeta;
    unsigned char umean;
    unsigned char vmean;
    char isom;
    int domx;
    int domy;
    list_node *next;
} node;

// Global variables
int N_BITALFA, N_BITBETA, min_size, max_size, SHIFT, image_width, image_height,
    int_max_alfa, isColor, MAX_ALFA, bits_per_coordinate_w,
    bits_per_coordinate_h, min, max, virtual_size;

FILE *input;
FILE *output;
node *trans;

// ----------------

long unpack(int size, FILE *fin) {
    int i;
    int value = 0;
    static int ptr = 1; /* how many bits are packed in sum so far */
    static int sum;

    /* size == -2 means we initialize things */
    if (size == -2) {
        sum = fgetc(fin);
        sum <<= 1;
        return ((long)0);
    }

    /* size == -1 means we want to peek at the next bit without */
    /* advancing the pointer */
    if (size == -1) return ((long)((sum & 256) >> 8));

    for (i = 0; i < size; ++i, ++ptr, sum <<= 1) {
        if (sum & 256) value |= 1 << i;

        if (ptr == 8) {
            sum = getc(fin);
            ptr = 0;
        }
    }
    return ((long)value);
}

int pack(int size, long value, FILE *foutf) {
    int i;
    static int ptr = 1, sum = 0, num_of_packed_bytes = 0;

    /* size == -1 means we are at the end, so write out what is left */
    if (size == -1 && ptr != 1) {
        fputc(sum << (8 - ptr), foutf);
        ++num_of_packed_bytes;
        return (0);
    }

    /* size == -2 means we want to know how many bytes we have written */
    if (size == -2) return (num_of_packed_bytes);

    for (i = 0; i<size; ++i, ++ptr, value = value>> 1, sum = sum << 1) {
        if (value & 1) sum |= 1;

        if (ptr == 8) {
            fputc(sum, foutf);
            ++num_of_packed_bytes;
            sum = 0;
            ptr = 0;
        }
    }
    return (-1);
}

void read_transformations(int atx, int aty, int size) {
    if (atx >= image_height || aty >= image_width) return;

    if (size > max_size || atx + size > image_height ||
        aty + size > image_width) {
        read_transformations(atx, aty, size / 2);
        read_transformations(atx + size / 2, aty, size / 2);
        read_transformations(atx, aty + size / 2, size / 2);
        read_transformations(atx + size / 2, aty + size / 2, size / 2);
        return;
    }

    if (size > min_size && unpack(1, input)) {
        /* A 1 means we subdivided.. so write 1 to output file and quadtree */
        // pack(1, (long)1, output);
        read_transformations(atx, aty, size / 2);
        read_transformations(atx + size / 2, aty, size / 2);
        read_transformations(atx, aty + size / 2, size / 2);
        read_transformations(atx + size / 2, aty + size / 2, size / 2);
    } else {
        /* Read the trasformation */
        // pack(1, (long)0, output);
        trans->next = (node *)malloc(sizeof(node));
        trans = trans->next;
        trans->next = NULL;
        trans->qalfa = (int)unpack(N_BITALFA, input);
        trans->qbeta = (int)unpack(N_BITBETA, input);
        if (isColor) {
            trans->umean = (int)unpack(8, input);
            trans->vmean = (int)unpack(8, input);
        }

        if (trans->qalfa != 0) {
            trans->isom = (int)unpack(3, input);
            trans->domx = SHIFT * (int)unpack(bits_per_coordinate_h, input);
            trans->domy = SHIFT * (int)unpack(bits_per_coordinate_w, input);
        } else {
            trans->isom = 0;
            trans->domx = 0;
            trans->domy = 0;
        }
        /* Write the domain coordinates to output */
        // pack(bits_per_coordinate_h, (long)trans->domx, output);
        // pack(bits_per_coordinate_w, (long)trans->domy, output);
    }
}

int main(int argc, char *argv[]) {
    char filein[100];
    char fileout[100];

    if (argc < 2) {
        printf("usage: ./pproc <input filename>\n");
        return -1;
    }

    // Open the input file
    strcpy(filein, argv[1]);
    FILE *input = fopen(filein, "r");

    if (input == NULL) {
        printf("error: could not read file %s\n", filein);
        return -1;
    }

    // Open the output file
    int i = 0;
    while (filein[i] != '.') {
        fileout[i] = filein[i];
        i++;
    }
    strcpy(&fileout[i], ".ifsp");

    FILE *output = fopen(fileout, "w");
    unpack(-2, input);

    // Read the header info first
    N_BITALFA = (int)unpack(4, input);
    N_BITBETA = (int)unpack(4, input);
    min_size = (int)unpack(7, input);
    max_size = (int)unpack(7, input);
    SHIFT = (int)unpack(6, input);
    image_width = (int)unpack(12, input);
    image_height = (int)unpack(12, input);
    int_max_alfa = (int)unpack(8, input);
    isColor = (int)unpack(1, input);

    // Read the transformations
    bits_per_coordinate_w = ceil(log(image_width / SHIFT) / log(2.0));
    bits_per_coordinate_h = ceil(log(image_height / SHIFT) / log(2.0));

    max = image_height;
    min = image_width;
    if (image_width > image_height) {
        min = image_height;
        max = image_width;
    }

    MAX_ALFA = (double)int_max_alfa / (double)(1 << 8) * (8.0);
    trans = (node *)malloc(sizeof(node));

    virtual_size = 1 << (int)ceil(log((double)max) / log(2.0));

    read_transformations(0, 0, virtual_size);

    // Close the files
    pack(-1, (long)0, output);
    pack(-2, (long)0, output);

    fclose(input);
    fclose(output);

    return 0;
}
