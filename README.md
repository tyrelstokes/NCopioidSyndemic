# NCopioidSyndemic
## Modeling of the opioid syndemic in North Carolina

In this project we estimate four latent factors characterizing variation of the syndemic across space and time in North Carolina. There are six different outcomes considered in this project: 
- Death counts: county-level counts; included are the counts of unintentional overdose deaths involving illicit opioids.
- ED visits: county-level counts; included are all ED visits related to medication and drug overdoses.
- Treatment: county-level counts; included are the number of uninsured individuals and Medicaid beneficiaries with opioid use disorder served by treatment programs.
- Buprenorphine: county-level counts; counts include residents who receive buprenorphine prescriptions.
- HCV: county-level counts; included acute and chronic cases of HCV cases.
- HIV: county-level counts; included are newly diagnosed HIV infections.

Description of the files:
- **NCdata.xlsx** includes the data;
- **Adj_Matrix_NC.RData** is the adjacency matrix used in our ICAR model;
- **NC_opioid.R** is the R-code we used to estimate the four factors.
