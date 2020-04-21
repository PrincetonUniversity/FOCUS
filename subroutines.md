# Subroutines list

There are produced pdf documentations for each individual file.
They contain some brief introductions.
Here are the list and links.

*Some of these might need to be updated. This is for the new version FOCUS. If you have questions, please contact Caoxiang Zhu.*

| names   | description    | source | documentation |
| -----   | -------------- | ------ | ------------- |
| globals | declaring global variables | [gloabls.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/globals.f90) | [globals.pdf](https://princetonuniversity.github.io/FOCUS/globals.pdf) |
| initial | initilization | [initial.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/initial.f90) | [initial.pdf](https://princetonuniversity.github.io/FOCUS/initial.pdf) |
| rdsurf  | read and discretize the unknot(toroidal) plasma boundary | [rdsurf.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/rdsurf.f90) | [rdsurf.pdf](https://princetonuniversity.github.io/FOCUS/rdsurf.pdf) |
| rdknot  | read and discretize the knotted boundary | [rdknot.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/rdknot.f90) | [rdknot.pdf](https://princetonuniversity.github.io/FOCUS/rdknot.pdf) |
| rdcoils | read and initilize the coils | [rdcoils.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/rdcoils.f90) | [rdcoils.pdf](https://princetonuniversity.github.io/FOCUS/rdcoils.pdf) |
| rdknot  | handling knotatron stuffs | [rdknot.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/rdknot.f90) | [rdknot.pdf](https://princetonuniversity.github.io/FOCUS/rdknot.pdf) |
| bfield  | Biot-Savart law for a current loop | [bfield.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/bfield.f90) | [bfield.pdf](https://princetonuniversity.github.io/FOCUS/bfield.pdf) |
| bnormal | evaluate the normal field residue | [bnormal.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/bnormal.f90) | [bnormal.pdf](https://princetonuniversity.github.io/FOCUS/bnormal.pdf) |
| bmnharm | evaluate the Bn spectrum difference | [bmnharm.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/bmnharm.f90) | [bmnharm.pdf](https://princetonuniversity.github.io/FOCUS/bmnharm.pdf) |
| torflux | evaluate the toroidal flux error | [torflux.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/torflux.f90) | [torflux.pdf](https://princetonuniversity.github.io/FOCUS/torflux.pdf) |
| length  | evaluate the coil length constraint | [length.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/length.f90) | [length.pdf](https://princetonuniversity.github.io/FOCUS/length.pdf) |
| surfsep | push coils away from any surfaces | [surfsep.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/surfsep.f90) | [surfsep.pdf](https://princetonuniversity.github.io/FOCUS/surfsep.pdf) |
| curvature| coil curvature metric|  [curvature.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/curvature.f90) | [curvature.pdf](https://princetonuniversity.github.io/FOCUS/curvature.pdf))
| solvers | interface for all the optimizers | [solvers.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/solvers.f90) | [solvers.pdf](https://princetonuniversity.github.io/FOCUS/solvers.pdf) |
| descent | differential flow optimization | [descent.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/descent.f90) | [descent.pdf](https://princetonuniversity.github.io/FOCUS/descent.pdf) |
| congrad | conjugate gradient optimization | [congrad.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/congrad.f90) | [congrad.pdf](https://princetonuniversity.github.io/FOCUS/congrad.pdf) |
| lmalg   | levenberg-marquardt optimization | [lmalg.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/lmalg.f90) | [lmalg.pdf](https://princetonuniversity.github.io/FOCUS/lmalg.pdf) |
| saving  | save the results | [saving.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/saving.f90) | [saving.pdf](https://princetonuniversity.github.io/FOCUS/saving.pdf) |
| focus   | main flow of the code | [focus.f90](https://github.com/PrincetonUniversity/FOCUS/blob/master/source/focus.f90) | [focus.pdf](https://princetonuniversity.github.io/FOCUS/focus.pdf) |
