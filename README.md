##Setup
Before running the program, install the required packages using the `requirements.txt` file. 

##Sequence extraction
To extract the gene sequences of the proteins, run the `preprocessing.py` script on the same directory containing the
`Supporting_Documents` folder.

```commandline
python preprocessing.py
```

##Robinson-Fould metric
To calculate the Robinson-Foulds metric, run the `robinson_foulds.py` script with the `Trees` folder containing in
the same directory.

```commandline
python robinson_foulds.py
```