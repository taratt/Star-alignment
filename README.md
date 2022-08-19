# Star-alignment
This project aims to implement the progressive multiple sequence alignment Star alignment for protein sequences with a block-based optimization and was done as an assignment for Introduction to Bioinformatics course at Amirkabir University of Technology.
## Scoring
The match score is simply considered to be equal to 3 and the mismatch and gap penalties are -1 and -2 respectively.
## Optimization
In this particular method, blocks containing more than or equal to two positions where each position has at least one mismatch or gap will be aligned and substituted by the original block one by one to see if the alignment score gets any higher. If the score rises by this action, the original block will be substituted by the align block.
## Input
The inputs to this code are 1) an integer which is the number of sequences that are going to be aligned together and 2) the protein sequences.
<br>
Example:
```
4 
TYIMREAQYESAQ
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES
```
## Output
The outputs to this code are 1) an integer which represents the score of the alignment and 2) The MSA.
<br>
Example:
```
48
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQE--R
WRYIAMRE-QYES--
```
