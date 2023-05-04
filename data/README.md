# Result data of Pathfinder
### **34 protein folding pathway with intermediate state structures**



**Data analysis:**   
                Zhaohong Huang
                College of Information Engineering  
                Zhejiang University of Technology, Hangzhou 310023, China  
		Email: 2112003269@zjut.edu.cn  

​                Xinyue Cui
​                College of Information Engineering  
​                Zhejiang University of Technology, Hangzhou 310023, China  
​		Email: 2112103122@zjut.edu.cn  

**Contact:**  
                Guijun Zhang, Prof  
                College of Information Engineering  
                Zhejiang University of Technology, Hangzhou 310023, China  
                Email: zgj@zjut.edu.cn  

## 1. File Directory

Pathfinder will eventually correct the seed state, not all seed states are intermediate states, so the seed state may not be included in the path. The number of intermediate states is adaptive, so different proteins have different intermediate states.

- data

  - Folding pathway
  - native.pdb
  - seedX.pdb 

- readme

- ContactOrder.py: Compute Contact Order at intermediate state and Contact Order at residue level. The residue-level Contact Order is appended to the b-factor of the pdb file by comparison with the native structure.

  ```
  $> python ContactOrder.py -p [file Directory]
  ```

  **note**: This script can calculate ContactOrder for all PDB files in this directory, but the native.pdb structure is required.

See the paper for details: Pathfinder: protein folding pathway prediction based on conformational sampling. Zhaohong Huang, Xinyue Cui, Yuhao Xia, Kailong Zhao, Guijun Zhang bioRxiv 2023.04.20.537604; doi [https://doi.org/10.1101/2023.04.20.537604](https://www.biorxiv.org/content/10.1101/2023.04.20.537604v1).




## 2. DISCLAIMER
The executable software and the source code of Pathfinder is distributed free of charge 
as it is to any non-commercial users. The authors hold no liabilities to the performance 
of the program.

