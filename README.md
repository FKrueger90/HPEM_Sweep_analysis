# HPEM_Sweep_analysis

This script is used to analyze parameter swweeps performed with the Hybrid Plasma Equipment Model (HPEM) develepod by Prof. Mark Kushner of the University of Michigan.

The script expects the directory structure to follow certain guidelines and wil not function properly if these guidelinies are violatged.

    .    
    ├── Base directory
    │   ├──  Case1 directory
    |   |    ├── icp.nam file
    |   |    ├── icp.log file
    |   |    ├── icp.out file 
    |   ├── Case2 directory
    |   |    ├── icp.nam file
    |   |    ├── icp.log file
    |   |    ├── icp.out file
    |   ├── ...   
    ├── Source directory (optional)
    ├── Template directory (optional) 

or 

    .    
    ├── Base directory
        ├── Case1 super directory
        |   ├── icp.nam file
        |   ├── ... 
        |   ├── Case1 directory
        |   |   ├── icp.nam file
        |   |   ├── icp.log file
        |   |   ├── icp.out file
        |   |   ├── ...
        ├── Case2 super directory
        |   ├── icp.nam file
        |   ├── ... 
        |   ├── Case2 directory
        |   |   ├── icp.nam file
        |   |   ├── icp.log file
        |   |   ├── icp.out file
        |   |   ├── ...
        ├── Source directory (optional)
        ├── Template directory (optional)
        ├── ...   

