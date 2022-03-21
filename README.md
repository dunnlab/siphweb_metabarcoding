# Characterizing the secret diets of siphonophores (Cnidaria: Hydrozoa) using DNA metabarcoding

Authors: Alejandro Damian-Serrano, Elizabeth D. Hetherington, C. Anela Choy, Steven H.D. Haddock, Alexandra Lapides, Casey W. Dunn

Corresponding author: Alejandro Damian-Serrano (email: adamians@uoregon.edu)

## Abstract:

Siphonophores (Cnidaria: Hydrozoa) are abundant and diverse gelatinous predators in open-ocean ecosystems. Due to limited access to the midwater, little is known about the diets of most deep-dwelling gelatinous species, which constrains our understanding of food-web structure and nutrient flow in these vast ecosystems. Visual gut-content methods can rarely identify soft-bodied rapidly-digested prey, while observations from submersibles often overlook small prey items. These methods have been differentially applied to shallow and deep siphonophore taxa, confounding habitat and methodological biases. DNA metabarcoding can be used to assess both shallow and deep species’ diets under a common methodological framework, since it can detect both small and gelatinous prey. We (1) further characterized the diets of open-ocean siphonophores using DNA metabarcoding, (2) compared the prey detected by visual and molecular methods to evaluate their technical biases, and (3) evaluated tentacle-based predictions of diet. To do this, we performed DNA metabarcoding analyses on the gut contents of 39 siphonophore species across depths to describe their diets, using six barcode regions along the 18S gene. Taxonomic identifications were assigned using public databases combined with local zooplankton sequences. We identified 55 unique prey items, including crustaceans, gelatinous animals, and fish across 47 siphonophore specimens in 24 species. We reported 29 novel predator-prey interactions, among them the first insights into the diets of nine siphonophore species, many of which were congruent with the dietary predictions based on tentilla morphology. Our analyses detected both small and gelatinous prey taxa underrepresented by visual methods in species from both shallow and deep habitats, indicating that siphonophores play similar trophic roles across depth habitats. We also reveal hidden links between siphonophores and filter-feeders near the base of the food web. This study expands our understanding of the ecological roles of siphonophores in the open ocean, their trophic roles within the ‘jelly-web’, and the importance of their diversity for nutrient flow and ecosystem functioning. Understanding these inconspicuous yet ubiquitous predator-prey interactions is critical to predict the impacts of climate change, overfishing, and conservation policies on oceanic ecosystems.

## Technical notes:


### Removed reference sequences:

To build one of our reference libraries, we collected a variety of potential prey items offshore from California, where most of the siphonophore samples were collected, using trawl nets and ROV samplers. From these specimens we extracted DNA, amplified the 18S+ITS-1 region using the two furthermost primer pairs used for metabarcoding, and sequenced them using Sanger. These were meant to supplement the scarce representation of deep midwater taxa in the SILVA databases. We submitted 89 sequences to GenBank in 2021, which became part of the SILVA138_Plus reference database we used to estimate taxonomic assignments from our Illumina gut content read data. 

In 2022, we ran a BLAST search on each of the sequences we submitted, and found that three of our chordate (two fishes and one salp) sequences are taxonomically mislabeled, as they correspond to calanoid copepod crustaceans. This is most likely due to the PCR amplifying DNA contamination from the environment or the tissue sample. 
The problematic accessions were:

MZ333597.1 Stenobrachius leucopsarus voucher 101 small subunit ribosomal RNA gene, partial sequence
MZ333623.1 Scopelogadus bispinosus voucher 284 small subunit ribosomal RNA gene, partial sequence
MZ333592.1 Salpa younti voucher 94 small subunit ribosomal RNA gene, partial sequence

These were removed from GenBank. Since these are only three sequences among millions in the reference database, we believe this misanotation would have no effect on the taxonomic assignment results. The fish, salp, and copepod prey sequences we found in our siphonophore gut contents were congruently identified among both reference databases used, with and without these reference sequences.

### Barcode region renaming:

We originally named the barcode regions by their predicted amplicon length (134, 152, 166, 179, 261, 272). For the PLOS ONE article, we renamed these  in the manuscript and figures to indicate the 18S regions they correspond to.

Here's the name conversion table:

| Original barcode names | GeneRegion  | Start position | End position | New barcode names |
| ---------------------- | ----------- | -------------- | ------------ | ----------------- |
| 134         | V9     |    1675   | 1790 | V9  |
| 152         | Between V5 and V7 (short amplicon)  |  1187  |  1339 |  V5-V7S  |
| 166         | V3          | 420 | 566 |   V3 |
| 179         | V7       |  1319   |  1489 | V7 |
| 261         | Part of V7 and all of V8      |  1472   |  1687  | V7p+V8  |
| 272         | Between V5 and V7 (long amplicon)  |  1067   | 1339 | V5-V7L  |
