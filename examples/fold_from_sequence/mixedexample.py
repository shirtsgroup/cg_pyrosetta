HETNAM     CG1 A   1  CG11                                                      
HETNAM     CG2 A   2  CG21                                                      
HETNAM     CG3 A   3  CG31                                                      
HETNAM     CG2 A   4  CG21                                                      
HETNAM     CG1 A   5  CG11                                                      
ATOM      1 BB1  CG1 A   1       0.000   0.000   0.000  1.00  0.00           X  
ATOM      2 SC1  CG1 A   1       0.773   0.131  -0.620  1.00  0.00           X  
ATOM      3 BB1  CG2 A   2      -0.500  -0.866  -0.000  1.00  0.00           X  
ATOM      4 BB2  CG2 A   2      -0.000  -1.732  -0.000  1.00  0.00           X  
ATOM      5 SC1  CG2 A   2      -1.500  -0.866   0.000  1.00  0.00           X  
ATOM      6 BB1  CG3 A   3       0.925  -1.775  -0.378  1.00  0.00           X  
ATOM      7 BB2  CG3 A   3       1.058  -1.852  -1.366  1.00  0.00           X  
ATOM      8 BB3  CG3 A   3       0.558  -1.249  -1.986  1.00  0.00           X  
ATOM      9 SC1  CG3 A   3       1.716  -1.742   0.233  1.00  0.00           X  
ATOM     10 BB1  CG2 A   4      -0.412  -1.064  -1.831  1.00  0.00           X  
ATOM     11 BB2  CG2 A   4      -1.049  -1.823  -1.690  1.00  0.00           X  
ATOM     12 SC1  CG2 A   4      -0.746  -0.122  -1.817  1.00  0.00           X  
ATOM     13 BB1  CG1 A   5      -0.709  -2.763  -1.666  1.00  0.00           X  
ATOM     14 SC1  CG1 A   5      -0.954  -3.401  -2.396  1.00  0.00           X  
TER                                                                             
CONECT    1    2    3                                                           
CONECT    2    1                                                                
CONECT    3    1    4    5                                                      
CONECT    4    3    6                                                           
CONECT    5    3                                                                
CONECT    6    4    7    9                                                      
CONECT    7    6    8                                                           
CONECT    8    7   10                                                           
CONECT    9    6                                                                
CONECT   10    8   11   12                                                      
CONECT   11   10   13                                                           
CONECT   12   10                                                                
CONECT   13   11   14                                                           
CONECT   14   13                                                                
# All scores below are weighted scores, not raw scores.
#BEGIN_POSE_ENERGIES_TABLE mixedexample.py
label fa_atr fa_rep fa_intra_atr fa_intra_rep mm_twist total
weights 1 1 1 1 0.1 NA
pose -4.43796 0.1047 -0.01642 0 0.12348 -4.22621
CG11_1 -1.06854 0.01216 0 0 0 -1.05638
CG21_2 -1.01572 0.04019 0 0 0.03737 -0.93815
CG31_3 -0.93384 0.00844 -0.01642 0 0.06174 -0.88007
CG21_4 -1.41987 0.0439 0 0 0.02436 -1.3516
CG11_5 0 0 0 0 0 0
#END_POSE_ENERGIES_TABLE mixedexample.py

