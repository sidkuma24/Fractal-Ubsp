Fractal Mars
============

       /\     Mars --  version 1.0 (10/28/1998)
      /__\    A quadtree based fractal image coder/decoder  
     /\  /\   Copyright (C) 1998 Mario Polvere <marpol@iname.com>     
    /__\/__\  University of Salerno Italy 

LICENSE
=======

This program is free software; you can redistribute it and/or 
modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation; either version 2 
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public 
License along with this program; if not, write to the Free 
Software Foundation, Inc., 59 Temple Place - Suite 330, 
Boston, MA 02111-1307, USA

GENERAL 
=======
This program is aimed to compare several speed-up 
techniques in fractal image coding. In particular six methods
have been implemented, chosen from both  classification  and 
feature vectors approaches. The methods are: Fisher, Hurtgen, 
MassCenter, Saupe, Saupe-Fisher and MassCenter-Saupe. Please
refer to the MANUAL for further details. 
I tried to keep it simple, many solutions have been adopted
for their simplicity and not because they are the best. 
Some fragments of code are repeted several times. This is because
I tried to make the various speed-up methods indipendent each
other in such a way you can see each method as a plug-in. So
if you want to try a new method you need just to write the code
that implement it and to insert it in the coder. In particular
you need to write only two functions.

Other features supported:
- Entropy based split decision function
- Variance based split decision function
- Adaptive splitting thresholds
- Output quadtree partition in pgm format
- Iterative decoding
- Piramidal decoding
- Fractal zooming


INSTALLATION
============

On most Unix systems you need just to type 
"make". You will get two binaries "encmars" and "decmars".

GETTING STARTED
===============

If you have a 512x512 raw image named
"lena.raw" in your working directory, you have just to
type "encmars" and the program will start with its default
configuration. You can also use an image n pgm format by 
typing "encmars image.pgm". Please to get a complete list 
of options type "encmars -h" and "decmars -h". 

KNOWN BUGS
==========

Strange lines may appear on the decoded image 
(on the block boundaries) when not power of two zooming 
factors are used.


AUTHOR
======

Mario Polvere <marpol@iname.com> University of
Salerno Italy. Comments, suggestions and bugs reports
are welcome.


ACKNOWLEDGEMENTS
================

I would like to thank Michele Nappi from
the University of Salerno who introduced me in the fractal
image coding world and made this work possible. Thanks also
to Yuval Fisher who wrote enc.c/dec.c from which some ideas
have been adopted and Matthias Ruhl who wrote frap/defrap 
from which kd-tree routines have been taken. 


