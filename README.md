# ComputAstroPhys_GroupReport
`main_with_openmp.cpp`: main code with PLM and (half)PPM. after running the code, you will get a output file `output.txt`. The file can be read and graphed by `test_C_result.py`.

## General Guidelines and Requirements
1. Must parallelize your program by at least one of the following methods \
a. OpenMP (minimum requirement) \
b. MPI (get bonus point) \
c. GPU (get extra bonus point)
2. Must provide convincing demonstration of the accuracy and performance of your program
3. Must use GitHub for collaborative development\
a. Tutorial: https://gitbook.tw/ \
b. Use the Fork→Pull→Push→Pull Request→Merge workflow\
c. Do NOT just upload the final code → must keep the development history on GitHub
4. Students per group: 3-4\
a. Inform the TA before April 26
5. Final presentation: May 31, June 7\
a. 25 mins presentation + 10 mins questions\
b. All students are encouraged to ask ANY questions
6. Bonus points: Surprise Me!!!
## Project-I: High-order Data Reconstruction in Hydro
1. Implement at least one of the following data reconstruction methods
a. Piecewise Parabolic Method (PPM)\
https://www.sciencedirect.com/science/article/pii/0021999184901438 \
https://www.sciencedirect.com/science/article/pii/S0021999108001435 \
b. Piecewise Cubic Method (PCM)\
https://www.sciencedirect.com/science/article/pii/S0021999117302759?via%3Dihub \
c. Weighted Essentially Non-Oscillatory (WENO)\
https://www.sciencedirect.com/science/article/pii/S0021999116301371 \
d. Gaussian Process (GP)\
https://link.springer.com/article/10.1007/s10915-017-0625-2 \
https://www.sciencedirect.com/science/article/pii/S0021999119300099 \
3. Compare the accuracy and performance with PLM
4. Bonus points \
a. Implement the characteristic tracing step (described in some of these schemes, e.g., PPM & PCM) to improve the temporal resolution\
b. Extend to 2D or even 3D \
c. MHD\
d. Add external and/or self-gravity
