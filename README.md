# Project Description
R Scripts and relating files. The aim of the project is to be able to process data collected from the single-molecule microscope and be able to detect anomalies within the collected data. The outcome of the project will be a decrease in the time taken to analyse data collected, hence increase efficiency and the ability to handle a larger quantity datasets. Main.r contains multiple functions that are key for analysis of the data.

# How to Use Scripts:
1. Before you intiallise the scripts ensure the following:
- Correct folder containing the datasets has been set to working directory e.g. use the command: setwd("LOCATION OF SETS")
- Ensure all required libraries are installed: ggplot, dplyr, pracma, gridExtra, and itsmr.
2. Inialise the script (press the source button or use the command source("")).
3. Use the command my.SMAgraph(175) to create a graph showing central tendency.

>Scripts
>- my.NormDat() - returns list of normalised dataframes
>- my.SMA() - returns list of smooth moving average dataframes
>- my.graph() - creates basic graphs
>- my.SMAgraph(k) - 
>- my.gradient(f) - 
>- my.Average() *Incomplete*


# Tasks



1. Load Csv - DONE
2. Normilase data by dividing the min y - DONE

3. Compare normalised data models e.g. statistical 

4. Plot norm graphs - DONE

5. Fit average Graphs if data is good
6. Plot average graphs


>TODO< 
>- Create a modle of the dat using function given by remi, start model at Xo (t=5).
>  Use the model to compare gradients, use the SMA.
>- Create a graph title function: use folder name as title e.g.
>  TacXPD-Y165E-35-Deg2-4.
