#Two distinct inhibitory neuronal classes govern acquisition and recall of spinal sensorimotor learning
##Figure replication

###Step 1: Data
All the data you need to run this cord are stored under:
> ConditioningExample/DataExp1/Day1/pair1

It is made of three raw kinematic files returning the position of 5 joints (Foot, ankle, knee, hip and crest) along the three axis X, Y and Z

###Step 2: Format the kinematic component of each joint and for Learner and Control
In order to separate Learner and Control kinematics and to remove the potential artefact from the dataset, run:
'''sh
sorting_raw_kinematic_file.m
'''
It will save the formatted kinematics in:
> ConditioningExample/DataExp1/AnalyseDay1/pair1.
