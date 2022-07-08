# Finite element solution of the generalized eigenvalue problem

An example code for the solution of the generalized eigenvalue problem in structural analysis. [Eigen library](https://eigen.tuxfamily.org/index.php?title=Main_Page) is used for linear algebra.

1. In the directory "example", an example of a vibration analysis given. The "problem_description.jpg" illustrates the problem considered, and the "mesh_slenderRing.txt" is the mesh used.
2. The incompatible mode element [1] is used for discretization, and the generalized eigenvalue problem is solved using Spectra library.
3. To run the executable "main" in the directory "build", open the terminal in the directory and enter

```
  ./main ../mesh_slenderRing.txt
```

## Reference
[1] Bathe, Klaus-Jürgen. Finite element procedures. Klaus-Jurgen Bathe, 2006.




