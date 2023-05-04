# Pathfinder
### **protein folding pathway prediction based on conformational sampling**



**Developer:**   
                Zhaohong Huang
                College of Information Engineering  
                Zhejiang University of Technology, Hangzhou 310023, China  
		Email: 2112003269@zjut.edu.cn  

**Contact:**  
                Guijun Zhang, Prof  
                College of Information Engineering  
                Zhejiang University of Technology, Hangzhou 310023, China  
                Email: zgj@zjut.edu.cn  

## 1. INSTALLATION
Binaries for Linux 64 bit system has been included in the package. The Linux binary was compiled using GCC 5.4.0. Users need to have these versions of GCC compilers when using binaries.

Please Follow the below steps to install and configure Pathfinder:

- Download Rosetta3.10 source package from https://www.rosettacommons.org/software/.
(where `$ROSETTA3_stage1`=path-to-Rosetta)

- Replicate the package (where `$ROSETTA3_stage2`=path-to-Rosetta_replocate), two Rosettas are used for stage 1 and stage 2 respectively.

- Copy and paste "tools" and "flags" from "code/" folder in Pathfinder package to work directory.

- Copy and paste ``"ClassicAbinitio.cc"`` and ``"ClassicAbinitio.hh"`` from ``"code/stage1/"`` folder in Pathfinder package to ``"$ROSETTA3_stage1/main/source/src/protocols/abinitio/"`` folder in Rosetta.

- Copy and paste ``"FragmentMover.cc"`` and ``"FragmentMover.hh"`` from ``"code/stage1/"`` folder in Pathfinder package to ``"$ROSETTA3_stage1/main/source/src/protocols/simple_moves/"`` folder in Rosetta.

- Copy and paste ``"ClassicAbinitio.cc"`` and ``"ClassicAbinitio.hh"`` from ``"code/stage2/"`` folder in Pathfinder package to ``"$ROSETTA3_stage2/main/source/src/protocols/abinitio/"`` folder in Rosetta.

- Copy and paste ``"FragmentMover.cc"`` and ``"FragmentMover.hh"`` from ``"code/stage2/"`` folder in Pathfinder package to ``"$ROSETTA3_stage2/main/source/src/protocols/simple_moves/"`` folder in Rosetta.
- Compile Pathfinder source code using the following commands:

```
 $> cd $ROSETTA3_stage1/main/source/
 $> ./scons.py AbinitioRelax -j<NumOfJobs> mode=release bin
 $> cd $ROSETTA3_stage2/main/source/
 $> ./scons.py AbinitioRelax -j<NumOfJobs> mode=release bin
```

## 2. INPUT
Pathfinder requires four files to generate models:

	fasta				: fasta file
	distance			: distance map file
	aat000_03_05.200_v1_3		: fragment library with fragment lenth 3
	aat000_09_05.200_v1_3		: fragment library with fragment lenth 9

## 3. RUN
Please follow the below steps to run Pathfinder:

- Go to the ``"example/"`` folder of Pathfinder.

- Run Pathfinder with the following command:

```
 $> $ROSETTA3_stage1/main/source/bin/AbinitioRelax.default.linuxgccrelease @flags 
 $> $ROSETTA3_stage2/main/source/bin/AbinitioRelax.default.linuxgccrelease @flags 
```

- Five models are generated in the ``"output_files/"`` folder.


## 4. OUTPUT
Output files of Pathfinder are stored in the ``"example/output_files/"`` folder, including self-adaptive predicted models (seed_X.pdb).

	data
	seed_1.pdb
	seed_2.pdb
	seed_3.pdb
	seed_4.pdb
	seed_5.pdb
	...


## 5. DISCLAIMER
The executable software and the source code of Pathfinder is distributed free of charge 
as it is to any non-commercial users. The authors hold no liabilities to the performance 
of the program.
