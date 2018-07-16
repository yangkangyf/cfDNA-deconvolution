## Deconvolution experiments 

* to make artifical mixture of cfDNA cpgs 
    * took ROADMAP differentially methylated WGBS data 
    * Selected some relevant tissues- brain, muscle- along with other tissues with a range of plausability- liver, heart atrium, cd34+ cells, 4star cell line etc 
    * Took random sample of 1000 of sites and the 1000 sites with greatest variance, decided to proceed with sites with greatest variance for now  
    * plotted the correlation of the methylation values - the mature organs seem to be the most similar (ie right atrium/right ventricle, thymus and spleen) and the most different seems to be the 4star cell line 
    * brain tissue not as similar as the two heart tissues
    * start by assuming 10 total reads 
* Mixture 1: average of two brain tissues for observed 
* Mixture 2: average of two brain tissues + random jitter @10% 
* Mixture 3: average of two brain tissues and thymus 

* TO DO: 

	* different mixture, different noise, different sites, tissue not in the reference, change %
	* quantify the error, do plots like before with simulation 
	* use the sites picked for mybaits with summing +/- 3 sites 
	* add noise to some of the tissues 
	* leave % of the sites missing, representative of the real data
	* what is better - adding +/-3 or 5 or 10 sites, versus just the only site   
	
	* real data: restrict to 450K sites, +/- summed RRBS/WGBS sites 
	
	* error for each site comes from read depth
	* imputation 
