HETNAM     CG1 A   1  CG11                                                      
HETNAM     CG1 A   2  CG11                                                      
HETNAM     CG1 A   3  CG11                                                      
HETNAM     CG1 A   4  CG11                                                      
HETNAM     CG1 A   5  CG11                                                      
ATOM      1 BB1  CG1 A   1       0.000   0.000   0.000  1.00  0.00           X  
ATOM      2 SC1  CG1 A   1      -0.439   0.831  -0.342  1.00  0.00           X  
ATOM      3 BB1  CG1 A   2      -0.500  -0.866  -0.000  1.00  0.00           X  
ATOM      4 SC1  CG1 A   2      -0.000  -1.732  -0.000  1.00  0.00           X  
ATOM      5 BB1  CG1 A   3      -1.500  -0.866  -0.000  1.00  0.00           X  
ATOM      6 SC1  CG1 A   3      -2.000  -1.480   0.611  1.00  0.00           X  
ATOM      7 BB1  CG1 A   4      -2.000  -0.253  -0.611  1.00  0.00           X  
ATOM      8 SC1  CG1 A   4      -2.788   0.259  -0.270  1.00  0.00           X  
ATOM      9 BB1  CG1 A   5      -1.712  -0.151  -1.563  1.00  0.00           X  
ATOM     10 SC1  CG1 A   5      -0.798  -0.446  -1.843  1.00  0.00           X  
TER                                                                             
CONECT    1    2    3                                                           
CONECT    2    1                                                                
CONECT    3    1    4    5                                                      
CONECT    4    3                                                                
CONECT    5    3    6    7                                                      
CONECT    6    5                                                                
CONECT    7    5    8    9                                                      
CONECT    8    7                                                                
CONECT    9    7   10                                                           
CONECT   10    9                                                                
# All scores below are weighted scores, not raw scores.
#BEGIN_POSE_ENERGIES_TABLE CG11example.py
label fa_atr fa_rep fa_intra_atr fa_intra_rep mm_twist total
weights 1 1 1 1 0.1 NA
pose -2.21948 0.40977 0 0 1.2 -0.60971
CG11_1 -0.63989 0.00792 0 0 0.3 -0.33196
CG11_2 -0.47001 0.19696 0 0 0.3 0.02696
CG11_3 0 0 0 0 0 0
CG11_4 -0.46532 0.19571 0 0 0.3 0.03039
CG11_5 -0.64427 0.00917 0 0 0.3 -0.3351
#END_POSE_ENERGIES_TABLE CG11example.py

