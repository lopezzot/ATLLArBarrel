# ATLLArBarrel
A Geant4 simulation of the ATLAS LAr Barrel beam test setup.

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#project-description">Project description</a></li>
    <li><a href="#authors-and-contacts">Authors and contacts</a></li>
    <li>
      <a href="#how-to">How to</a>
      <ul>
        <li><a href="#build-compile-and-execute-on-maclinux">Build, compile and execute on Mac/Linux</a></li>
        <li><a href="#build-compile-and-execute-on-lxplus">Build, compile and execute on lxplus</a></li>
      </ul>
    </li>
  </ol>
</details>

<!--Project desription-->
## Project description
The project targets a standalone Geant4 simulation of the ATLAS LAr Barrel test beam setup. It is used for Geant4 regression testing and physics lists comparison, as well as for testing speeding up solutions.
- ⏰ Start date: February 2023
- 📌 Status: development

<!--Authors and contacts-->
## Authors and contacts
- 👨‍🔬 Lorenzo Pezzotti (CERN EP-SFT) - lorenzo.pezzotti@cern.ch 
- 👨‍🔬 Supervisor: Alberto Ribon (CERN EP-SFT)

<!--How to-->
## How to

### Build, compile and execute on Mac/Linux
1.  git clone the repo
    ```sh
    git clone https://github.com/lopezzot/ATLLArBarrel.git
    ```
2.  source Geant4-11.1 env
    ```sh
    source /relative_path_to/geant4.11.1-install/bin/geant4.sh
    ```
3.  cmake build directory and make (using geant4.11.1)
    ```sh
    mkdir build; cd build/
    cmake -DGeant4_DIR=/absolute_path_to/geant4.11.1-install/lib/Geant4-11.1/ relative_path_to/ATLLarBarrel/
    make
    ```
    See CMake options for all build options
4.  execute (example with run.mac macro card, 2 threads and FTFP_BERT physics list)
    ```sh
    ./ATLLarBarrel -m run.mac -t 2 -p FTFP_BERT
    ```

Parser options
- `-m macro.mac`: pass a Geant4 macro card (example `-m run.mac` available in source directory and automatically copied in build directory) 
- `-t integer`: pass number of threads for multi-thread execution (example `-t 2`, default is the number of threads on the machine)
- `-p Physics_List`: select Geant4 physics list (example `-p FTFP_BERT`)

### Build, compile and execute on lxplus
1. git clone the repo
   ```sh
   git clone https://github.com/lopezzot/ATLLArBarrel.git
   ```
2. cmake build directory and make (using geant4.11.1, check for gcc and cmake dependencies for other versions)
   ```sh
   mkdir build; cd build/
   source /cvmfs/sft.cern.ch/lcg/contrib/gcc/8.3.0/x86_64-centos7/setup.sh
   source /cvmfs/geant4.cern.ch/geant4/11.1/x86_64-centos7-gcc8-optdeb-MT/CMake-setup.sh
   export CXX=`which g++`
   export CC=`which gcc`
   cmake3 -DGeant4_DIR= /cvmfs/geant4.cern.ch/geant4/11.1/x86_64-centos7-gcc8-optdeb-MT/lib64/Geant4-11.1.0/ ../ATLLArBarrel/
   ```
3. execute (example with TBrun.mac macro card, 2 threads and FTFP_BERT physics list)
   ```sh
   ./ATLLArBarrel -m run.mac -t 2 -p FTFP_BERT
   ```
