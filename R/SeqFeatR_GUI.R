#R GUI written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-07-27 to 2012-?-?
#GUI for scripts to better put data in it

#search for "change" for points where changes are possible
#please check before executing
rm(list = ls())

#libs
require(tcltk)
require(tcltk2)
require(widgetTools)
#----------------------------------------------------------------------------------------------------------------
#loaded objects
calculate_pos_epi_is_running <- FALSE
calculate_pos_epi_further_is_running <- FALSE
calculate_co_mut_is_running <- FALSE
calculate_q_value_is_running <- FALSE
calculate_co_mut_graphics_is_running <- FALSE
calculate_co_mut_graphics_both_is_running <- FALSE
calculate_founder_is_running <- FALSE
create_sequence_graphic_is_running <- FALSE
shared_mutations_is_running <- FALSE
rewrite_shared_mutations_is_running <- FALSE

#----------------------------------------------------------------------------------------------------------------
#Error handling---------------------------------------------------------------------------------------------------
Require <- structure(function(
	### checks if a package is available
	pkg
	### the package to check for
	){
    	if (data.class(result<-try(find.package(pkg),TRUE))=="try-error")    {
        	tkmessageBox(title="An error has occured!",message=paste("Cannot find package\"",pkg,"\". Please install it with install.packages(\"",pkg,"\")",sep=""),icon="error",type="ok")
        	return (FALSE)
    	}else{
        	require(pkg,character.only=TRUE)
       		return (TRUE)
    	}
	### True if package is there, False if it is not
},ex=function(){
	Require(tcltk)
})


#Helper-Functions for the program start---------------------------------------------------------------------------
get_label_right <- function(splitted){
	length <- length(splitted[[1]])
	if (length>1){
		return (paste(splitted[[1]][2], "/../", splitted[[1]][length]))
	}
	else{
		return (splitted[[1]][2])
	}
}

getcsvfile <- function(labelT, funct, whichfile)  {
	loaded_files <- .GlobalEnv[["loaded_files"]]
     	name <- tclvalue(tkgetOpenFile(filetypes="{{csv Files} {.csv}} {{All files} *}"))
	loaded_files[1, whichfile, funct] <- name
	splitted <- strsplit(name, "/")
	label <- get_label_right(splitted)
	tclvalue(labelT) <- label
	.GlobalEnv[["loaded_files"]] <- loaded_files
}

getFASTAfile <- function(labelT, funct, whichfile)  {
	loaded_files <- .GlobalEnv[["loaded_files"]]
     	name <- tclvalue(tkgetOpenFile(filetypes="{{FASTA Files} {.fasta}} {{All files} *}"))
	loaded_files[1, whichfile, funct] <- name
	splitted <- strsplit(name, "/")
	label <- get_label_right(splitted)
	tclvalue(labelT) <- label
	.GlobalEnv[["loaded_files"]] <- loaded_files
}

getnexusfile <- function(labelT, funct, whichfile)  {
	loaded_files <- .GlobalEnv[["loaded_files"]]
     	name <- tclvalue(tkgetOpenFile(filetypes="{{nexus Files} {.nh}} {{All files} *}"))
	loaded_files[1, whichfile, funct] <- name
	splitted <- strsplit(name, "/")
	label <- get_label_right(splitted)
	tclvalue(labelT) <- label
	.GlobalEnv[["loaded_files"]] <- loaded_files
}

check_phylo <- function(cbValue, textEntryWidgetKIB){
	if (cbValue=="1"){
		tkconfigure(textEntryWidgetKIB, state="normal")
	}else {
		tkconfigure(textEntryWidgetKIB, state="disabled")
	}
}

check_allels <- function(cbValue, cMtextEntryWidgetC, button.widget_co2, textEntryWidgetHC, textEntryWidgetH1C, textEntryWidgetH2C, textEntryWidgetH12C, textEntryWidgetHBC, textEntryWidgetH1BC, textEntryWidgetH2BC, textEntryWidgetH12BC) {
	if (cbValue=="1"){
		tkconfigure(cMtextEntryWidgetC, state="normal")
		tkconfigure(button.widget_co2, state="disabled")
		tkconfigure(textEntryWidgetHC, state="normal")
		tkconfigure(textEntryWidgetH1C, state="normal")
		tkconfigure(textEntryWidgetH2C, state="normal")
		tkconfigure(textEntryWidgetH12C, state="normal")
		tkconfigure(textEntryWidgetHBC, state="normal")
		tkconfigure(textEntryWidgetH1BC, state="normal")
		tkconfigure(textEntryWidgetH2BC, state="normal")
		tkconfigure(textEntryWidgetH12BC, state="normal")
	}else {
		tkconfigure(cMtextEntryWidgetC, state="disabled")
		tkconfigure(button.widget_co2, state="normal")
		tkconfigure(textEntryWidgetHC, state="disabled")
		tkconfigure(textEntryWidgetH1C, state="disabled")
		tkconfigure(textEntryWidgetH2C, state="disabled")
		tkconfigure(textEntryWidgetH12C, state="disabled")
		tkconfigure(textEntryWidgetHBC, state="disabled")
		tkconfigure(textEntryWidgetH1BC, state="disabled")
		tkconfigure(textEntryWidgetH2BC, state="disabled")
		tkconfigure(textEntryWidgetH12BC, state="disabled")
	}
}

check_dna <- function(cbValue) {
	if (cbValue=="1"){
		
	}else {
		
	}
}

check_ident <- function(cbValue, textEntryWidgetH, textEntryWidgetH1, textEntryWidgetH12, textEntryWidgetH2, textEntryWidgetHB, textEntryWidgetH1B, textEntryWidgetH12B, textEntryWidgetH2B, textEntryWidgetHC, textEntryWidgetH1C, textEntryWidgetH2C, textEntryWidgetH12C, textEntryWidgetHBC, textEntryWidgetH1BC, textEntryWidgetH2BC, textEntryWidgetH12BC, smtextEntryWidgetX, smtextEntryWidgetH, smtextEntryWidgetH1, smtextEntryWidgetH2, smtextEntryWidgetH12, smtextEntryWidgetHB, smtextEntryWidgetH1B, smtextEntryWidgetH2B, smtextEntryWidgetH12B) {
	if (cbValue=="1"){
		tkconfigure(textEntryWidgetH, state="disabled")
		tkconfigure(textEntryWidgetH1, state="disabled")
		tkconfigure(textEntryWidgetH12, state="disabled")
		tkconfigure(textEntryWidgetH2, state="disabled")
		tkconfigure(textEntryWidgetHB, state="disabled")
		tkconfigure(textEntryWidgetH1B, state="disabled")
		tkconfigure(textEntryWidgetH12B, state="disabled")
		tkconfigure(textEntryWidgetH2B, state="disabled")
		tkconfigure(textEntryWidgetHC, state="disabled")
		tkconfigure(textEntryWidgetH1C, state="disabled")
		tkconfigure(textEntryWidgetH2C, state="disabled")
		tkconfigure(textEntryWidgetH12C, state="disabled")
		tkconfigure(textEntryWidgetHBC, state="disabled")
		tkconfigure(textEntryWidgetH1BC, state="disabled")
		tkconfigure(textEntryWidgetH2BC, state="disabled")
		tkconfigure(textEntryWidgetH12BC, state="disabled")
		tkconfigure(smtextEntryWidgetX, state="normal")
		tkconfigure(smtextEntryWidgetH, state="disabled")
		tkconfigure(smtextEntryWidgetH1, state="disabled")
		tkconfigure(smtextEntryWidgetH2, state="disabled")
		tkconfigure(smtextEntryWidgetH12, state="disabled")
		tkconfigure(smtextEntryWidgetHB, state="disabled")
		tkconfigure(smtextEntryWidgetH1B, state="disabled")
		tkconfigure(smtextEntryWidgetH2B, state="disabled")
		tkconfigure(smtextEntryWidgetH12B, state="disabled")
	}else {
		tkconfigure(textEntryWidgetH, state="normal")
		tkconfigure(textEntryWidgetH1, state="normal")
		tkconfigure(textEntryWidgetH12, state="normal")
		tkconfigure(textEntryWidgetH2, state="normal")
		tkconfigure(textEntryWidgetHB, state="normal")
		tkconfigure(textEntryWidgetH1B, state="normal")
		tkconfigure(textEntryWidgetH12B, state="normal")
		tkconfigure(textEntryWidgetH2B, state="normal")
		tkconfigure(textEntryWidgetHC, state="normal")
		tkconfigure(textEntryWidgetH1C, state="normal")
		tkconfigure(textEntryWidgetH2C, state="normal")
		tkconfigure(textEntryWidgetH12C, state="normal")
		tkconfigure(textEntryWidgetHBC, state="normal")
		tkconfigure(textEntryWidgetH1BC, state="normal")
		tkconfigure(textEntryWidgetH2BC, state="normal")
		tkconfigure(textEntryWidgetH12BC, state="normal")
		tkconfigure(smtextEntryWidgetX, state="disabled")
		tkconfigure(smtextEntryWidgetH, state="normal")
 		tkconfigure(smtextEntryWidgetH1, state="normal")
 		tkconfigure(smtextEntryWidgetH2, state="normal")
 		tkconfigure(smtextEntryWidgetH12, state="normal")
 		tkconfigure(smtextEntryWidgetHB, state="normal")
 		tkconfigure(smtextEntryWidgetH1B, state="normal")
 		tkconfigure(smtextEntryWidgetH2B, state="normal")
 		tkconfigure(smtextEntryWidgetH12B, state="normal")
	}
}

generate_example_fasta <- function(){
	sequences <- c("LPDIQGNENMGYQPSWIFCGMETNGSQCLEEMFHCCWINC", "MPDWNQKWGNDHLASINLD-WLKTIQQPGIEKHLRFYENW", "VPDASGKHGIIGMDVTSSMERRHGMVQLPWPAMVWGRPHW", "MPDVRGVGCARRDCLIVHRFCMPFNNQVYCKVWIVYWTYK", "QPDTPKITRKEATAIHKCGIHWQTNCQKLSTVHPFHHQVD", "SWDDFSDFTMVHQWYAQGTLGPYKAMQLKMIFQGVSIMEV", "IPDEPCYCCVKNKILTVEIGVHHAKSQVRRNIDNIRRKTE", "HFST-ICPYIWKMYFTWMGQKLVIQKVNGRTPPHCDECNQ", "SNFT-TTKLRDQHNLYPAGLQEIEHKVDHQILGIYGQIWY", "ETSTALRTQDQTFMLALRANYMVMLKVLDCISVKLFICWR", "DSSTMDAECSTLQRFIWWHAHYAWIRVAKKPYCLDCPYAV", "KKSTLGIARGIQRSHGWYWRQTHCVMVLTPSQHKMGEKSW", "ICSTELCGCLINWPPMQWIVFAHMDDVNDSQTNTCDMRSQ", "GPSTNARTMGGQDCAYMTHTLTKHIWVILAFDPIMIVHKP")
	names <- c("P01_HLA_A01_00_B01_02", "P02_HLA_A01_00_B01_02", "P03_HLA_A01_02_B01_02", "P04_HLA_A01_00_B01_02", "P05_HLA_A01_00_B01_02", "P06_HLA_A01_02_B01_02", "P07_HLA_A01_02_B01_02", "P08_HLA_A04_03_B04_03", "P09_HLA_A04_03_B04_03", "P10_HLA_A04_03_B04_03", "P11_HLA_A04_03_B04_00", "P12_HLA_A04_03_B04_03", "P13_HLA_A04_03_B04_00", "P14_HLA_A04_03_B04_03")
	aa_seqs <- AAStringSet(sequences)
	names(aa_seqs) <- names
	save_file <- tclvalue(tkgetSaveFile(initialfile = "Example_aa.fasta", filetypes = "{{Fasta Files} {.fasta}} {{All files} *}"))
	writeXStringSet(aa_seqs, save_file)
}

generate_example_config <- function(){
	config <- matrix(c("Number of patients threshold", "Height of horizontal bar", "Height of star level", "Level for possible epitope", "Window Threshold", "HLA-A1", "HLA-A2", "HLA-A3", "HLA-A4", "HLA-B1", "HLA-B2", "HLA-B3", "HLA-B4", "p-value correction", "Phylogenetic comparison", "Matrix for phylo", "One identifier", "window_size", "DNA", "p-value addon", "csv seperator", "number of cases", "position of odds.ratio", "position of p.values", "name for the y-axis label", "frequency", "position of amino acid", "maximum value of y-axis", "intervall of y-axis", "add color information", "Number of patients threshold_cm", "level for significance", "More than one core", "cmHLA-A1", "cmHLA-A2", "cmHLA-A3", "cmHLA-A4", "cmHLA-B1", "cmHLA-B2", "cmHLA-B3", "cmHLA-B4", "cmp-value correction", "Ratio in subtree", "cMlevel for significance", "space between blocks", "colors of the plot", "name positions", "names", "ticks", "column of first position in first file", "column of second position in first file", "column of values in first file", "column of first position in second file", "column of second position in second file", "column of values in second file", "threshold", "min_number_of_ele_in_tupel", "max_number_of_ele_in_tupel", "column", "column of values", "column of aas", "smHLA-A1", "smHLA-A2", "smHLA-A3", "smHLA-A4", "smHLA-B1", "smHLA-B2", "smHLA-B3", "smHLA-B4", "Identifier", "With Allels", "cbWith Allels", "rw column of first position", "rw column of second position", "rw column of values", "rw csv seperator", "rw threshold", "1", "0.01", "0.5", "0.01", "0.01", "10", "11", "13", "14", "17", "18", "20", "21", "bonferroni", "0", "BLOSUM62", "0", "9", "0", "0.01", ";", "8", "6", "2", "2", "7", "NULL", "4", "0.2", "0", "1", "0.01", "0", "10", "11", "13", "14", "17", "18", "20", "21", "bonferroni", "7", "0.01", "5", "wheat, darkblue, black, green", "1,3", "S, F", "1,5", "4", "5", "7", "2", "3", "4", "0.2", "2", "2", "9", "9", "12", "10", "11", "13", "14", "17", "18", "20", "21", "X4", "0", "0", "3", "4", "10", "\\t", "0.1"), ncol=2)
	print (config)
	save_file <- tclvalue(tkgetSaveFile(initialfile = "example_config.cfg", filetypes = "{{Config Files} {.cfg}} {{All files} *}"))
	write.table(config, save_file, append = FALSE, quote = FALSE, sep = "=", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")
}

save_config <- function(z){
	table <- c()
	save_file <- tclvalue(tkgetSaveFile(initialfile = "user_config.cfg", filetypes = "{{Config Files} {.cfg}} {{All files} *}"))
	for (i in 1:length(z)){
		key <- z[[i]][[2]]
		value <- tclvalue(z[[i]][[1]])
		table <- rbind(table, c(key, value))
	}
	write.table(table, save_file, append = FALSE, quote = FALSE, sep = "=", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")
}

read.config <- function(config){
	tryCatch({
		key.val <- read.table(config, sep="=", col.names=c("key","value"), as.is=c(1,2))
		return (key.val)
	}, error = function(err) {tkmessageBox(title="Error",message=paste("No config file loaded. Standard parameters will be used!\n", err),icon="error",type="ok")})
	
}

attach.config <- function(config){
	z <- .GlobalEnv[["z"]]
	tryCatch({
		con <- read.table(config, sep="=", col.names=c("key","value"), as.is=c(1,2))
		for (i in 1:length(z)){
			object <- z[[i]][[1]]
			object_name <- z[[i]][[2]]
			nopt <- con[which(con[,1]==object_name),2]
			tclvalue(object) <- nopt
		}
	}, error = function(err) {tkmessageBox(title="Error",message=paste("Error in config file.\n", err),icon="error",type="ok")})
}

open_tutorial <- function(){
	epi <- vignette("SeqFeatR")
	if (length(epi) > 0){
		print (epi)
	}
}

about <- function(){
	tryCatch({
		ReturnVal <- tkmessageBox(title = "About SeqFeatR", message = paste("SeqFeatR Version ", packageVersion("SeqFeatR"), sep=""), icon = "info", type = "ok")
	}, error = function(err) {tkmessageBox(title="Error",message=paste("SeqFeatR is not installed"),icon="error",type="ok")})
}

#Functions for the program start-------------------------------------------------------------------------------------------
calculate_pos_epi <- function(tbn, cbValue01, cbValue_0e, textEntryC, textEntryD, textEntryE, textEntryF, textEntryG, textEntryH, textEntryH1, textEntryH2, textEntryH12, textEntryHB, textEntryH1B, textEntryH2B, textEntryH12B, textEntryIB, textEntryKIB, txt, icbValue01, textEntryH42B){
	Require("Biostrings")
	Require("plyr")
	Require("plotrix")
	loaded_files <- .GlobalEnv[["loaded_files"]]
	tkconfigure(tbn,cursor="watch")
	if (!calculate_pos_epi_is_running){
		result <- c()
		calculate_pos_epi_is_running <- TRUE
		if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
			source("assocpoint.R")
		}
		tryCatch({ep_wo_set_input_file_known_sequences(loaded_files[1,1,1])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No sequence file",icon="error",type="ok")})
		save_name_epi <- tclvalue(tkgetSaveFile(initialfile="epitope_results.pdf",filetypes="{{PDF Files} {.pdf}} {{All files} *}"))
			if (!nchar(save_name_epi)){
    				tkmessageBox(message="No file was selected!")
			}
		save_name_epi_csv <- tclvalue(tkgetSaveFile(initialfile="epitope_results.csv",filetypes="{{csv Files} {.csv}} {{All files} *}"))
			if (!nchar(save_name_epi_csv)){
    				tkmessageBox(message="No file was selected!")
			}
		if(as.numeric(tclvalue(textEntryH42B)) <=2){
			print ("Beware! Window size smaller than 2! Will be set to 2")
			tclvalue(textEntryH42B) <- 2
		}
		tryCatch({result <- assocpoint(loaded_files[1,1,1], loaded_files[1,2,1], loaded_files[1,3,1], save_name_epi, save_name_epi_csv, as.numeric(tclvalue(cbValue01)), as.numeric(tclvalue(textEntryC)), as.numeric(tclvalue(textEntryD)), as.numeric(tclvalue(textEntryE)), as.numeric(tclvalue(textEntryF)), as.numeric(tclvalue(textEntryG)), as.numeric(tclvalue(textEntryH)), as.numeric(tclvalue(textEntryH1)), as.numeric(tclvalue(textEntryH2)), as.numeric(tclvalue(textEntryH12)), as.numeric(tclvalue(textEntryHB)), as.numeric(tclvalue(textEntryH1B)), as.numeric(tclvalue(textEntryH2B)), as.numeric(tclvalue(textEntryH12B)), tclvalue(textEntryIB), as.numeric(tclvalue(cbValue_0e)), tclvalue(textEntryKIB), loaded_files[1,4,1], as.logical(as.numeric(icbValue01)), as.numeric(tclvalue(textEntryH42B)) )}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		old.option <- (getOption("width", default = NULL))
		options(width=2000)
		text <- paste(capture.output(print(result)),collapse="\n")
		options(width=old.option)
		tkconfigure(txt, state="normal")
		tryCatch({tkinsert(txt,"end",text)}, error = function(err) {})
		tkconfigure(txt, state="disabled")
		calculate_pos_epi_is_running <- FALSE
	}
	tkconfigure(tbn,cursor="arrow")
}

calculate_pos_epi_further <- function(tbn, textEntry_p.value, txt){
	Require("Biostrings")
	loaded_files <- .GlobalEnv[["loaded_files"]]
	tkconfigure(tbn,cursor="watch")
	if (!calculate_pos_epi_further_is_running){
		result_ef <- c()
		calculate_pos_epi_further_is_running <- TRUE
		if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
			source("assocpointtuple.R")
		}
    		tryCatch({pos_epi_get_input_file_sequences(loaded_files[1,1,4])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No sequence file",icon="error",type="ok")})
		tryCatch({pos_epi_get_input_file_pos_epi(loaded_files[1,2,4])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No possible epitope file",icon="error",type="ok")})
		tryCatch({pos_epi_get_input_file_mut_core(loaded_files[1,3,4], as.numeric(tclvalue(textEntry_p.value)))}, error = function(err) {tkmessageBox(title="An error has occured!",message="No mutation core file",icon="error",type="ok")})
		save_name_ef1 <- tclvalue(tkgetSaveFile(initialfile="co_mut_in_epitopes_result.csv",filetypes="{{csv Files} {.csv}} {{All files} *}"))
			if (!nchar(save_name_ef1)){
    				tkmessageBox(message="No file was selected!")
			}
		save_name_ef2 <- tclvalue(tkgetSaveFile(initialfile="possible_compensatory_mutation.csv",filetypes="{{csv Files} {.csv}} {{All files} *}"))
			if (!nchar(save_name_ef2)){
    				tkmessageBox(message="No file was selected!")
			}
		tryCatch({result_ef <- assocpointtuple(loaded_files[1,1,4], loaded_files[1,2,4], loaded_files[1,3,4], as.numeric(tclvalue(textEntry_p.value)), save_name_ef1, save_name_ef2)}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		old.option <- (getOption("width", default = NULL))
		options(width=2000)
		text <- paste(capture.output(print (result_ef)),collapse="\n")
		options(width=old.option)
		tkconfigure(txt, state="normal")
		tryCatch({tkinsert(txt,"end",text)}, error = function(err) {})
		tkconfigure(txt, state="disabled")
		calculate_pos_epi_further_is_running <- FALSE
	}
	tkconfigure(tbn,cursor="arrow")
}

create_sequence_graphic_G <- function(tbn, BtextEntryA, BtextEntryB, BtextEntryC, BtextEntryD, BtextEntryE, BtextEntryF, BtextEntryG, BtextEntryH, BtextEntryI, txt, BtextEntryJ){
	Require("plotrix")
	tkconfigure(tbn,cursor="watch")
	loaded_files <- .GlobalEnv[["loaded_files"]]
	if (as.character(tclvalue(BtextEntryA)) == "\\t"){
		sep <- "\t"
	}
	else{
		sep <- as.character(tclvalue(BtextEntryA))
	}
	if (!create_sequence_graphic_is_running){
		result <- c()
		create_sequence_graphic_is_running <- TRUE
		if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
			source("orPlot.R")
		}
		tryCatch({two_wo_set_input_file_known_sequences(loaded_files[1,1,7], sep)}, error = function(err) {tkmessageBox(title="An error has occured!",message="No result input file",icon="error",type="ok")})
		save_name_two <- tclvalue(tkgetSaveFile(initialfile="final_graphic",filetypes="{{PDF Files} {.pdf}} {{All files} *}"))
		if (!nchar(save_name_two)){
    			tkmessageBox(message="No file was selected!")
		}
		tryCatch({result <- orPlot(loaded_files[1,1,7], save_name_two, sep, as.numeric(tclvalue(BtextEntryB)), as.numeric(tclvalue(BtextEntryC)), as.numeric(tclvalue(BtextEntryD)), as.numeric(tclvalue(BtextEntryE)), as.numeric(tclvalue(BtextEntryF)), as.numeric(tclvalue(BtextEntryG)), as.numeric(tclvalue(BtextEntryH)), as.numeric(tclvalue(BtextEntryI)),  as.numeric(tclvalue(BtextEntryJ)))}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		old.option <- (getOption("width", default = NULL))
		options(width=2000)
		text <- paste(capture.output(print(result)),collapse="\n")
		options(width=old.option)
		tkconfigure(txt, state="normal")
		print (text)
		tryCatch({tkinsert(txt,"end",text)}, error = function(err) {})
		tkconfigure(txt, state="disabled")
		create_sequence_graphic_is_running <- FALSE
	}
	tkconfigure(tbn,cursor="arrow")
}

calculate_co_mut <- function(allels, tbn2, txt, cMtextEntryC, cMtextEntryD, cmtextEntryIB, textEntryHC, textEntryH1C, textEntryH2C, textEntryH12C, textEntryHBC, textEntryH1BC, textEntryH2BC, textEntryH12BC, textEntryIB, cbValue01){
	Require("Biostrings")
	loaded_files <- .GlobalEnv[["loaded_files"]]
	tkconfigure(tbn2,cursor="watch")
	if (!calculate_co_mut_is_running){
		result_co <- c()
		calculate_co_mut_is_running <- TRUE
		if (allels=="1"){
			if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
				source("assocpairfeat.R")
			}
			tryCatch({cm_get_input_file(loaded_files[1,1,2])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No sequence file",icon="error",type="ok")})
			save_name_co_mut <- tclvalue(tkgetSaveFile(initialfile="co_mutation_results.csv",filetypes="{{csv Files} {.csv}} {{All files} *}"))
			if (!nchar(save_name_co_mut)){
    				tkmessageBox(message="No file was selected!")
			}
			tryCatch({result_co <- assocpairfeat(loaded_files[1,1,2], save_name_co_mut, as.numeric(tclvalue(cbValue01)), as.numeric(tclvalue(cMtextEntryC)), as.numeric(tclvalue(cMtextEntryD)), as.numeric(tclvalue(textEntryHC)), as.numeric(tclvalue(textEntryH1C)), as.numeric(tclvalue(textEntryH2C)), as.numeric(tclvalue(textEntryH12C)), as.numeric(tclvalue(textEntryHBC)), as.numeric(tclvalue(textEntryH1BC)), as.numeric(tclvalue(textEntryH2BC)), as.numeric(tclvalue(textEntryH12BC)), tclvalue(cmtextEntryIB))}, error = function(err) {tkmessageBox(title="An error has occured!",message="No comutation analysis done",icon="error",type="ok")})
		}	
		if (allels=="0"){
			if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
				source("assocpair.R")
			}
			tryCatch({cm_wo_get_input_file_sequences(loaded_files[1,1,2])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No sequence file",icon="error",type="ok")})
			tryCatch({cm_wo_get_input_file_consensus(loaded_files[1,2,2])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No consensus file",icon="error",type="ok")})
			save_name_co_mut <- tclvalue(tkgetSaveFile(initialfile="co_mutation_results_wo_allels.csv",filetypes="{{csv Files} {.csv}} {{All files} *}"))
			if (!nchar(save_name_co_mut)){
    				tkmessageBox(message="No file was selected!")
			}
			tryCatch({result_co <- assocpair(loaded_files[1,1,2], loaded_files[1,2,2], save_name_co_mut, as.numeric(tclvalue(cbValue01)), as.numeric(tclvalue(cMtextEntryD)), tclvalue(cmtextEntryIB))}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		}
		old.option <- (getOption("width", default = NULL))
		options(width=2000)
		text <- paste(capture.output(print (result_co)),collapse="\n")
		options(width=old.option)
		tkconfigure(txt, state="normal")
		tryCatch({tkinsert(txt,"end",text)}, error = function(err) {})
		tkconfigure(txt, state="disabled")
		calculate_co_mut_is_running <- FALSE
	}
	tkconfigure(tbn2,cursor="arrow")
}

calculate_co_mut_graphics <- function(allels, tbn2, cMtextEntryE){
	Require("Biostrings")
	loaded_files <- .GlobalEnv[["loaded_files"]]
	tkconfigure(tbn2,cursor="watch")
	if (!calculate_co_mut_graphics_is_running){
		calculate_co_mut_graphics_is_running <- TRUE
		if (allels=="1"){
			if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
				source("vispairfeat.R")
			}
			tryCatch({co_g_get_input_file(loaded_files[1,1,5])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})
			save_name  <- tclvalue(tkgetSaveFile(initialfile="co_mut_g_results.pdf",filetypes="{{PDF Files} {.pdf}} {{All files} *}"))
			if (!nchar(save_name)){
    				tkmessageBox(message="No file was selected!")
			} 
			tryCatch({vispairfeat(loaded_files[1,1,5], save_name, as.numeric(tclvalue(cMtextEntryE)))}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		}	
		if (allels=="0"){
			if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
				source("vispair.R")
			}
			tryCatch({co_g_get_input_file_wo_allel(loaded_files[1,1,5])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})
			save_name <- tclvalue(tkgetSaveFile(initialfile="co_mut_g_results_wo_allels.pdf",filetypes="{{PDF Files} {.pdf}} {{All files} *}"))
			if (!nchar(save_name)){
    				tkmessageBox(message="No file was selected!")
			}
			tryCatch({vispair(loaded_files[1,1,5], save_name, as.numeric(tclvalue(cMtextEntryE)))}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		}
		calculate_co_mut_graphics_is_running <- FALSE
	}
	tkconfigure(tbn2,cursor="arrow")
}

calculate_co_mut_graphics_both <- function(tbn5, txt, TcMtextEntryE, TcMtextEntryE1, TcMtextEntryE2, TcMtextEntryE3, TcMtextEntryE4, TcMtextEntryE5, TcMtextEntryE6, TcMtextEntryE7, TcMtextEntryE8, TcMtextEntryE9, TcMtextEntryE0){
	Require("Biostrings")
	loaded_files <- .GlobalEnv[["loaded_files"]]
	tkconfigure(tbn5,cursor="watch")
	if (!calculate_co_mut_graphics_both_is_running){
		calculate_co_mut_graphics_both_is_running <- TRUE
		if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
			source("tartan.R")
		}
		if (loaded_files[1,3,9] == "0"){
			with_distance_matrix <- FALSE
		} else{
			with_distance_matrix <- TRUE
		}
		tryCatch({co_g_get_input_file_wo_allel(loaded_files[1,1,9])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})
		tryCatch({co_g_get_input_file2_wo_allel(loaded_files[1,2,9])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})
		tryCatch({co_g_get_distance_matrix(loaded_files[1,3,9], with_distance_matrix)}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})
		save_name  <- tclvalue(tkgetSaveFile(initialfile="co_mut_g_both_results.pdf",filetypes="{{PDF Files} {.pdf}} {{All files} *}"))
		if (!nchar(save_name)){
    			tkmessageBox(message="No file was selected!")
		}
		list1 <- as.numeric(strsplit(as.character(tclvalue(TcMtextEntryE2)), ",")[[1]])
		list2 <- as.numeric(strsplit(as.character(tclvalue(TcMtextEntryE4)), ",")[[1]])
		tryCatch({result_tartan <- tartan(loaded_files[1,1,9], loaded_files[1,2,9], save_name, as.numeric(tclvalue(TcMtextEntryE)), strsplit(as.character(tclvalue(TcMtextEntryE1)), ",")[[1]], list1, strsplit(as.character(tclvalue(TcMtextEntryE3)), ",")[[1]], list2, as.numeric(tclvalue(TcMtextEntryE5)), as.numeric(tclvalue(TcMtextEntryE6)), as.numeric(tclvalue(TcMtextEntryE7)), as.numeric(tclvalue(TcMtextEntryE8)), as.numeric(tclvalue(TcMtextEntryE9)), as.numeric(tclvalue(TcMtextEntryE0)), with_distance_matrix, loaded_files[1,3,9] )}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		old.option <- (getOption("width", default = NULL))
		options(width=2000)
		text <- paste(capture.output(print (result_tartan)),collapse="\n")
		options(width=old.option)
		tkconfigure(txt, state="normal")
		tryCatch({tkinsert(txt,"end",text)}, error = function(err) {})
		tkconfigure(txt, state="disabled")
		calculate_co_mut_graphics_both_is_running <- FALSE
	}
	tkconfigure(tbn5,cursor="arrow")
}

calculate_founder <- function(tbn2, fetextEntryE, txt){
	Require("Biostrings")
	Require("phangorn")
	loaded_files <- .GlobalEnv[["loaded_files"]]
	tkconfigure(tbn2,cursor="watch")
	if (!calculate_founder_is_running){
		result_fe <- c()
		if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
			source("founderelim.R")
		}
		calculate_founder_is_running <- TRUE
		tryCatch({fe_wo_get_input_file_known_sequences(loaded_files[1,1,6])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})
		tryCatch({fe_wo_get_input_file_co_mut_results(loaded_files[1,2,6])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})
		tryCatch({fe_wo_get_input_file_tree(loaded_files[1,3,6])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})		
		save_name_fe <- tclvalue(tkgetSaveFile(initialfile="co_mut_results_from_founder.csv",filetypes="{{csv Files} {.csv}} {{All files} *}"))
			if (!nchar(save_name_fe)){
    				tkmessageBox(message="No file was selected!")
			}
    		tryCatch({result_fe <- founderelim(loaded_files[1,1,6], loaded_files[1,2,6], loaded_files[1,3,6], save_name_fe, as.numeric(tclvalue(fetextEntryE)))}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		old.option <- (getOption("width", default = NULL))
		options(width=2000)
		text <- paste(capture.output(print (result_fe)),collapse="\n")
		options(width=old.option)
		tkconfigure(txt, state="normal")
		tryCatch({tkinsert(txt,"end",text)}, error = function(err) {})
		tkconfigure(txt, state="disabled")
		calculate_founder_is_running <- FALSE
	}
	tkconfigure(tbn2,cursor="arrow")
}

calculate_q_value <- function(tbn3, txt){
	Require("Biostrings")
	Require("qvalue")
	loaded_files <- .GlobalEnv[["loaded_files"]]
	tkconfigure(tbn3,cursor="watch")
	if (!calculate_q_value_is_running){
		result_q <- c()
		if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
			source("q-values.R")
		}
		calculate_q_value_is_running <- TRUE
		tryCatch({q_value_set_input_file(loaded_files[1,1,3])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})		
		save_name_q <- tclvalue(tkgetSaveFile(initialfile="csv_with_q_values.csv",filetypes="{{csv Files} {.csv}} {{All files} *}"))
			if (!nchar(save_name_q)){
    				tkmessageBox(message="No file was selected!")
			}
    		tryCatch({result_q <- qvalues(loaded_files[1,1,3], save_name_q)}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		old.option <- (getOption("width", default = NULL))
		options(width=2000)
		text <- paste(capture.output(print (result_q, width="200")),collapse="\n")
		options(width=old.option)
		tkconfigure(txt, state="normal")
		tryCatch({tkinsert(txt,"end",text)}, error = function(err) {})
		tkconfigure(txt, state="disabled")
		calculate_q_value_is_running <- FALSE
	}
	tkconfigure(tbn3,cursor="arrow")
}

shared_mutations <- function(tbn5, txt, smtextEntryA, smtextEntryB, smtextEntryC, smtextEntryD, smtextEntryE, smtextEntryF, smtextEntryH, smtextEntryH1, smtextEntryH2, smtextEntryH12, smtextEntryHB, smtextEntryH1B, smtextEntryH2B, smtextEntryH12B, icbValue01, smtextEntryidet){
	Require("Biostrings")
	Require("parallel")
	loaded_files <- .GlobalEnv[["loaded_files"]]
	tkconfigure(tbn5,cursor="watch")
	if (!calculate_q_value_is_running){
		result_q <- c()
		if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
			source("assoctuple.R")
		}
		shared_mutations_is_running <- TRUE
		tryCatch({sm_wo_set_input_file_ori_sequences(loaded_files[1,1,8])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})
		tryCatch({sm_wo_set_input_file_epi_results(loaded_files[1,2,8])}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})	
		save_name_sm <- tclvalue(tkgetSaveFile(initialfile="shared_mutations.csv",filetypes="{{csv Files} {.csv}} {{All files} *}"))
		if (!nchar(save_name_sm)){
    			tkmessageBox(message="No file was selected!")
		}
    		tryCatch({result_sm <- assoctuple(loaded_files[1,1,8], loaded_files[1,2,8], as.numeric(tclvalue(smtextEntryA)), as.numeric(tclvalue(smtextEntryB)), as.numeric(tclvalue(smtextEntryC)), save_name_sm, as.numeric(tclvalue(smtextEntryD)), as.numeric(tclvalue(smtextEntryE)), as.numeric(tclvalue(smtextEntryF)), as.numeric(tclvalue(smtextEntryH)), as.numeric(tclvalue(smtextEntryH1)), as.numeric(tclvalue(smtextEntryH2)), as.numeric(tclvalue(smtextEntryH12)), as.numeric(tclvalue(smtextEntryHB)), as.numeric(tclvalue(smtextEntryH1B)), as.numeric(tclvalue(smtextEntryH2B)), as.numeric(tclvalue(smtextEntryH12B)), as.logical(as.numeric(icbValue01)), as.character(tclvalue(smtextEntryidet)) )}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		old.option <- (getOption("width", default = NULL))
		options(width=2000)
		text <- paste(capture.output(print (read.csv2(result_sm), width="200")),collapse="\n")
		options(width=old.option)
		tkconfigure(txt, state="normal")
		tryCatch({tkinsert(txt,"end",text)}, error = function(err) {})
		tkconfigure(txt, state="disabled")
		shared_mutations_is_running <- FALSE
	}
	tkconfigure(tbn5,cursor="arrow")
}

rewrite_shared_mutations <- function(tbn3, rsmtextEntryE5, rsmtextEntryE6, rsmtextEntryE7, rsmtextEntryA, rsmtextEntryB, txt){
	loaded_files <- .GlobalEnv[["loaded_files"]]
	tkconfigure(tbn3,cursor="watch")
	if (as.character(tclvalue(rsmtextEntryA)) == "\\t"){
		sep <- "\t"
	}
	else{
		sep <- as.character(tclvalue(rsmtextEntryA))
	}
	if (!rewrite_shared_mutations_is_running){
		result_q <- c()
		if (data.class(result<-try(find.package("SeqFeatR"),TRUE))=="try-error"){
			source("rewritetuple.R")
		}
		rewrite_shared_mutations_is_running <- TRUE
		tryCatch({rsm_wo_set_input_file_epi_results(loaded_files[1,1,10],  sep)}, error = function(err) {tkmessageBox(title="An error has occured!",message="No input file",icon="error",type="ok")})
		save_name_rsm <- tclvalue(tkgetSaveFile(initialfile="rewritten_shared_mutations.csv",filetypes="{{csv Files} {.csv}} {{All files} *}"))
		if (!nchar(save_name_rsm)){
    			tkmessageBox(message="No file was selected!")
		}
    		tryCatch({result_rsm <- rewrite_shared_mutations_result_inner(loaded_files[1,1,10], save_name_rsm, as.numeric(tclvalue(rsmtextEntryE5)), as.numeric(tclvalue(rsmtextEntryE6)), as.numeric(tclvalue(rsmtextEntryE7)), sep, as.numeric(tclvalue(rsmtextEntryB)))}, error = function(err) {tkmessageBox(title="An error has occured!",message=paste("No calculation done\n", err),icon="error",type="ok")})
		old.option <- (getOption("width", default = NULL))
		options(width=2000)
		text <- paste(capture.output(print(result_rsm)),collapse="\n")
		options(width=old.option)
		tkconfigure(txt, state="normal")
		tryCatch({tkinsert(txt,"end",text)}, error = function(err) {})
		tkconfigure(txt, state="disabled")
		rewrite_shared_mutations_is_running <- FALSE
	}
	tkconfigure(tbn3,cursor="arrow")
}

#----------------------------------------------------------
SeqFeatR_GUI <- structure(function(
	### SeqFeatR GUI to handle the program without command line
	){
loaded_files <- array(rep(0, 40), dim=c(1,4,10), dimnames = list("name",c("specific file_1","specific file_2","specific file_3", "specific file_4"),c("pos_epi", "co-mut", "q-value", "pos_e_cal", "graphics_co_mut", "founder_elim", "two_side", "shared_mutations", "co_mut_graphics_both", "rewrite_shared_mutations")))
.GlobalEnv[["loaded_files"]] <- loaded_files

if(length(system.file("extdata", "config.cfg", package="SeqFeatR")) > 1){
	con <- read.config(system.file("extdata", "config.cfg", package="SeqFeatR"))
}else {
	con <- c()
}

tt <- tktoplevel()

topMenu <- tkmenu(tt)           # Create a menu
tkconfigure(tt, menu = topMenu)
fileMenu <- tkmenu(topMenu, tearoff = FALSE)
helpMenu <- tkmenu(topMenu, tearoff = FALSE)

openRecentMenu <- tkmenu(topMenu, tearoff = FALSE)
tkadd(fileMenu, "command", label = "Open config file", command = function() attach.config(tclvalue(tkgetOpenFile(initialfile = "config.cfg", filetypes = "{{Config Files} {.cfg}} {{All files} *}"))))
tkadd(fileMenu, "command", label = "Save config file", command = function() save_config(z))
tkadd(fileMenu, "command", label = "Generate example fasta", command = function() generate_example_fasta())
tkadd(fileMenu, "command", label = "Generate example config", command = function() generate_example_config())
tkadd(fileMenu, "command", label = "Quit", command = function() tkdestroy(tt))

tkadd(helpMenu, "command", label = "Tutorial", command = function() open_tutorial())
tkadd(helpMenu, "command", label = "About SeqFeatR", command = function() about())

tkadd(topMenu, "cascade", label = "File", menu = fileMenu)
tkadd(topMenu, "cascade", label = "Help", menu = helpMenu)

frameOverall <- tkframe(tt)

frame_progs <- tkframe(frameOverall,relief="groove",borderwidth=2, background="white")
frame_info <- tkframe(frameOverall,relief="groove",borderwidth=2)

tkpack(frameOverall, fill="both",expand="yes")

tkgrid(frame_progs, sticky="snew")
tkgrid(frame_info)
#------------------------------------------------------------------------------------------------------------------------------
#tn <- tkwidget(frame_progs, "ttk::notebook")

#print (.Tcl.args(background="red"))

tcl("ttk::style", "configure", "TNotebook", background="white")
# Make non selected tabs a little darker
tcl("ttk::style", "configure", "TNotebook.Tab", background="lightgrey")
#tcl("ttk::style", "map", "TNotebook.Tab", background=c("active", "skyblue"))
#tcl("ttk::style", "map", "TNotebook.Tab", background=c("active", "red"))
tcl("ttk::style", "map", "TNotebook.Tab", background=c("selected","khaki1"))

tn <- ttknotebook(frame_progs)

#tkwm.resizable(tt, TRUE, TRUE)### FOR WINDOWS COMMENT OUT!!
tkwm.title(tt, "SeqFeatR - Epitope and Co-mutation")
tkgrid.columnconfigure(tt,0,weight=1)
tkgrid.rowconfigure(tt,0,weight=1)
tkwm.minsize(tt, "640", "480")

tkgrid.rowconfigure(frameOverall,1,weight=1, minsize=100)
tkgrid.rowconfigure(frameOverall,0,weight=10, minsize=200)
tkgrid.columnconfigure(frameOverall,0,weight=1, minsize=100)

#tkgrid.configure(frameOverall,sticky='nswe')
tkgrid.rowconfigure(frame_progs,0,weight=1, minsize=20)
tkgrid.rowconfigure(frame_progs,1,weight=100, minsize=100)
tkgrid.columnconfigure(frame_progs,0,weight=1, minsize=10)
tkgrid.configure(frame_progs,sticky='nswe')
tkgrid.rowconfigure(frame_info,1,weight=1, minsize=10)
tkgrid.columnconfigure(frame_info,1,weight=1, minsize=10)

#-----------------------------------------------------------------------------------------------------------------------------
dna <- tkframe(frame_progs)
tkgrid(dna)
tkgrid.rowconfigure(dna,0,weight=0)
tkgrid.columnconfigure(dna,0,weight=0)

cb01 <- tkcheckbutton(dna,command=function()check_dna(cbValu01 <- as.character(tclvalue(cbValue01))), background="#FF6600")
cbValue01 <- tclVar("0")
tkconfigure(cb01,variable=cbValue01)

icb01 <- tkcheckbutton(dna,command=function()check_ident(icbValu01 <- as.character(tclvalue(icbValue01)), textEntryWidgetH, textEntryWidgetH1, textEntryWidgetH12, textEntryWidgetH2, textEntryWidgetHB, textEntryWidgetH1B, textEntryWidgetH12B, textEntryWidgetH2B, textEntryWidgetHC, textEntryWidgetH1C, textEntryWidgetH2C, textEntryWidgetH12C, textEntryWidgetHBC, textEntryWidgetH1BC, textEntryWidgetH2BC, textEntryWidgetH12BC, smtextEntryWidgetX, smtextEntryWidgetH, smtextEntryWidgetH1, smtextEntryWidgetH2, smtextEntryWidgetH12, smtextEntryWidgetHB, smtextEntryWidgetH1B, smtextEntryWidgetH2B, smtextEntryWidgetH12B), background="#FF6600")
icbValue01 <- tclVar("0")
tkconfigure(icb01,variable=icbValue01)

tkgrid(tklabel(dna,text="Fasta files with nucleotides?", background="#FF6600"), cb01, tklabel(dna,text="Only one feature?", background="#FF6600"), icb01 )

#-----------------------------------------------------------------------------------------------------------------------------

tkgrid.rowconfigure(tn,0,weight=5)
tkgrid.columnconfigure(tn,0,weight=5)
tkgrid.configure(tn,sticky='nswe')

otbn <- ttkframe(tn)
otbn2 <- ttkframe(tn)
otbn3 <- ttkframe(tn)
otbn5 <- ttkframe(tn)

tkadd(tn,otbn,text="Calculate point mutations vs. feature(s)")   ### tabid=0
tkadd(tn,otbn2,text="Calculate tuples vs. feature(s)")  ### tabid=1
tkadd(tn,otbn5,text="Create graphics for tuples calculations") ### tabid=5
tkadd(tn,otbn3,text="Add-on scripts") ### tabid=3

#------------------------------------------------------------------------------------------------------------------------------
scr <- tkscrollbar(frame_info, repeatinterval=5, command=function(...)tkyview(txt,...))
xscr <- tkscrollbar(frame_info, repeatinterval=5, command=function(...)tkxview(txt,...),orient='horiz')
txt <- tktext(frame_info,bg="white", font="courier", height=15, width=80, wrap="none", yscrollcommand=function(...)tkset(scr,...), xscrollcommand=function(...)tkset(xscr,...))
tkgrid(txt,scr)
tkgrid(xscr)

tkgrid.configure(scr,sticky="ns")
tkgrid.configure(xscr,sticky="ew")

tkconfigure(txt, state="disabled")
tkgrid.rowconfigure(txt,0,weight=1)
tkgrid.columnconfigure(txt,0,weight=1)
tkgrid.configure(txt,sticky='nswe')

#---------------------------------------------------
tkgrid.rowconfigure(otbn,0,weight=1, minsize=100)
tkgrid.columnconfigure(otbn,0,weight=10, minsize=200)
tkgrid.columnconfigure(otbn,1,weight=3, minsize=100)

tbn <- ttkframe(otbn)
helptbn <- ttkframe(otbn)
tkgrid(tbn, helptbn, sticky="snew")
tkgrid.columnconfigure(tbn,0,weight=10)
tkgrid.columnconfigure(helptbn,0,weight=1)

scr1 <- tkscrollbar(tbn, command=function(...)tkyview(txt1,...),bg='white')
txt1 <- tktext(tbn,height=1, width=33, bg='grey') 
xscr1 <- tkscrollbar(tbn, command=function(...)tkxview(txt1,...),orient='horiz')
tkconfigure(txt1, yscrollcommand=function(...)tkset(scr1,...), xscrollcommand=function(...)tkset(xscr1,...), state="disabled")
tkpack(scr1,side='right',fill='y')
tkpack(txt1,expand='yes',fill='both')
tkpack(xscr1, side="bottom", fill="x")
sf1 <- tkcanvas(txt1, background="red")
tkwindow.create(txt1, 'end', window=sf1)

epi <- tkframe(sf1,relief="groove",borderwidth=2)
epist <- tkframe(sf1,relief="groove",borderwidth=2)

tkgrid(epi, sticky="snew")
tkgrid(epist, sticky="snew")

lb1 <- tklabel(epi, text=paste("Discover associations between point mutations and features"), background="#FFCC99")
tkgrid (lb1, sticky="snew", columnspan=8)

button.widget_epi1 <- tkbutton(epi,text="Select FASTA file",command=function()getFASTAfile(pElabelText1, 1, 1))
pElabelText <- tclVar("Selected sequences to analyse: ")
pElabel1 <- tklabel(epi,text=tclvalue(pElabelText))
tkconfigure(pElabel1,textvariable=pElabelText)

pElabelText1 <- tclVar("")
pElabel2 <- tklabel(epi,text=tclvalue(pElabelText1))
tkconfigure(pElabel2,textvariable=pElabelText1)
tkgrid(pElabel1, pElabel2, button.widget_epi1)

button.widget_epi2 <- tkbutton(epi,text="Select csv file",command=function()getcsvfile(pElabelText1A, 1, 2))
pElabelTextA <- tclVar("Selected known epitopes: ")
pElabel1A <- tklabel(epi,text=tclvalue(pElabelTextA))
tkconfigure(pElabel1A,textvariable=pElabelTextA)

pElabelText1A <- tclVar("")
pElabel2A <- tklabel(epi,text=tclvalue(pElabelText1A))
tkconfigure(pElabel2A,textvariable=pElabelText1A)
tkgrid(pElabel1A, pElabel2A, button.widget_epi2)

button.widget_epi3 <- tkbutton(epi,text="Select csv file",command=function()getcsvfile(pElabelText1B, 1, 3))
pElabelTextB <- tclVar("Selected known binding motifs: ")
pElabel1B <- tklabel(epi,text=tclvalue(pElabelTextB))
tkconfigure(pElabel1B,textvariable=pElabelTextB)

button.widget_epi4 <- tkbutton(epi,text="Select FASTA File",command=function()getFASTAfile(pElabelText1D, 1, 4))
pElabelTextD <- tclVar("Selected reference sequence: ")
pElabel1D <- tklabel(epi,text=tclvalue(pElabelTextD))
tkconfigure(pElabel1D,textvariable=pElabelTextD)

pElabelText1D <- tclVar("")
pElabel2D <- tklabel(epi,text=tclvalue(pElabelText1D))
tkconfigure(pElabel2D,textvariable=pElabelText1D)
tkgrid(pElabel1D, pElabel2D, button.widget_epi4)

pElabelText1B <- tclVar("")
pElabel2B <- tklabel(epi,text=tclvalue(pElabelText1B))
tkconfigure(pElabel2B,textvariable=pElabelText1B)
tkgrid(pElabel1B, pElabel2B, button.widget_epi3)

pElabelTextC <- tclVar("Minimal number of members: ")
pElabel1C <- tklabel(epi,text=tclvalue(pElabelTextC))
if (length(con[which(con[,1]=="Number of patients threshold"),2]) > 0){
	textEntryC <- tclVar(con[which(con[,1]=="Number of patients threshold"),2])
}else{
	textEntryC <- tclVar("")
}
textEntryWidgetC <- tkentry(epi,width=5,textvariable=textEntryC)
tkgrid(pElabel1C, textEntryWidgetC)

pElabelTextD <- tclVar("Height of horizontal bar: ")
pElabel1D <- tklabel(epi,text=tclvalue(pElabelTextD))
if (length(con[which(con[,1]=="Height of horizontal bar"),2]) > 0){
	textEntryD <- tclVar(con[which(con[,1]=="Height of horizontal bar"),2])
}else{
	textEntryD <- tclVar("")
}
textEntryWidgetD <- tkentry(epi,width=5,textvariable=textEntryD)
tkgrid(pElabel1D, textEntryWidgetD)

pElabelTextE <- tclVar("Height of star level: ")
pElabel1E <- tklabel(epi,text=tclvalue(pElabelTextE))
if (length(con[which(con[,1]=="Height of star level"),2]) > 0){
	textEntryE <- tclVar(con[which(con[,1]=="Height of star level"),2])
}else{
	textEntryE <- tclVar("")
}
textEntryWidgetE <- tkentry(epi,width=5,textvariable=textEntryE)
tkgrid(pElabel1E, textEntryWidgetE)

pElabelTextF <- tclVar("Level for possible association: ")
pElabel1F <- tklabel(epi,text=tclvalue(pElabelTextF))
if (length(con[which(con[,1]=="Level for possible epitope"),2]) > 0){
	textEntryF <- tclVar(con[which(con[,1]=="Level for possible epitope"),2])
}else{
	textEntryF <- tclVar("")
}
textEntryWidgetF <- tkentry(epi,width=5,textvariable=textEntryF)
tkgrid(pElabel1F, textEntryWidgetF)

pElabelTextG <- tclVar("Window threshold: ")
pElabel1G <- tklabel(epi,text=tclvalue(pElabelTextG))
if (length(con[which(con[,1]=="Window Threshold"),2]) > 0){
	textEntryG <- tclVar(con[which(con[,1]=="Window Threshold"),2])
}else{
	textEntryG <- tclVar("")
}
textEntryWidgetG <- tkentry(epi,width=5,textvariable=textEntryG)
tkgrid(pElabel1G, textEntryWidgetG)

pElabelTextH <- tclVar("Position of feature A: ")
pElabel1H <- tklabel(epi,text=tclvalue(pElabelTextH))
if (length(con[which(con[,1]=="HLA-A1"),2]) > 0){
	textEntryH <- tclVar(con[which(con[,1]=="HLA-A1"),2])
}else{
	textEntryH <- tclVar("")
}
textEntryWidgetH <- tkentry(epi,width=3,textvariable=textEntryH)
pElabelTextH1 <- tclVar("-")
pElabel1H1 <- tklabel(epi,text=tclvalue(pElabelTextH1))
if (length(con[which(con[,1]=="HLA-A2"),2]) > 0){
	textEntryH1 <- tclVar(con[which(con[,1]=="HLA-A2"),2])
}else{
	textEntryH1 <- tclVar("")
}
textEntryWidgetH1 <- tkentry(epi,width=3,textvariable=textEntryH1)
if (length(con[which(con[,1]=="HLA-A3"),2]) > 0){
	textEntryH2 <- tclVar(con[which(con[,1]=="HLA-A3"),2])
}else{
	textEntryH2 <- tclVar("")
}
textEntryWidgetH2 <- tkentry(epi,width=3,textvariable=textEntryH2)
pElabelTextH12 <- tclVar("-")
pElabel1H12 <- tklabel(epi,text=tclvalue(pElabelTextH12))
if (length(con[which(con[,1]=="HLA-A4"),2]) > 0){
	textEntryH12 <- tclVar(con[which(con[,1]=="HLA-A4"),2])
}else{
	textEntryH12 <- tclVar("")
}
textEntryWidgetH12 <- tkentry(epi,width=3,textvariable=textEntryH12)
tkgrid(pElabel1H, textEntryWidgetH, pElabel1H1, textEntryWidgetH1, textEntryWidgetH2, pElabel1H12, textEntryWidgetH12)

pElabelTextHB <- tclVar("Position of feature B: ")
pElabel1HB <- tklabel(epi,text=tclvalue(pElabelTextHB))
if (length(con[which(con[,1]=="HLA-B1"),2]) > 0){
	textEntryHB <- tclVar(con[which(con[,1]=="HLA-B1"),2])
}else{
	textEntryHB <- tclVar("")
}
textEntryWidgetHB <- tkentry(epi,width=3,textvariable=textEntryHB)
pElabelTextH1B <- tclVar("-")
pElabel1H1B <- tklabel(epi,text=tclvalue(pElabelTextH1B))
if (length(con[which(con[,1]=="HLA-B2"),2]) > 0){
	textEntryH1B <- tclVar(con[which(con[,1]=="HLA-B2"),2])
}else{
	textEntryH1B <- tclVar("")
}
textEntryWidgetH1B <- tkentry(epi,width=3,textvariable=textEntryH1B)
if (length(con[which(con[,1]=="HLA-B3"),2]) > 0){
	textEntryH2B <- tclVar(con[which(con[,1]=="HLA-B3"),2])
}else{
	textEntryH2B <- tclVar("")
}
textEntryWidgetH2B <- tkentry(epi,width=3,textvariable=textEntryH2B)
pElabelTextH12B <- tclVar("-")
pElabel1H12B <- tklabel(epi,text=tclvalue(pElabelTextH12B))
if (length(con[which(con[,1]=="HLA-B4"),2]) > 0){
	textEntryH12B <- tclVar(con[which(con[,1]=="HLA-B4"),2])
}else{
	textEntryH12B <- tclVar("")
}
textEntryWidgetH12B <- tkentry(epi,width=3,textvariable=textEntryH12B)
tkgrid(pElabel1HB, textEntryWidgetHB, pElabel1H1B, textEntryWidgetH1B, textEntryWidgetH2B, pElabel1H12B, textEntryWidgetH12B)

pElabelTextIB <- tclVar("P-value correction: ")
pElabel1IB <- tklabel(epi,text=tclvalue(pElabelTextIB))
ddb <- tkframe(epi)
textEntryIB <- tclVar()
dropdownList(ddb, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), textEntryIB, 15, "none")
tkgrid(pElabel1IB, ddb)

cb0e <- tkcheckbutton(epi,command=function()check_phylo(cbValu_0e <- as.character(tclvalue(cbValue_0e)), textEntryWidgetKIB), background="#FF6600")
if (length(con[which(con[,1]=="Phylogenetic comparison"),2]) > 0){
	cbValue_0e <- tclVar(con[which(con[,1]=="Phylogenetic comparison"),2])
}else{
	cbValue_0e <- tclVar("0")
}
tkconfigure(cb0e,variable=cbValue_0e)
tkgrid(tklabel(epi,text="Add phylogenetic analysis?", background="#FF6600"),cb0e, sticky="snew", columnspan=4)

pElabelTextKIB <- tclVar("Matrix for phylogenetic analysis: ")
pElabel1KIB <- tklabel(epi,text=tclvalue(pElabelTextKIB))
if (length(con[which(con[,1]=="Matrix for phylo"),2]) > 0){
	textEntryKIB <- tclVar(con[which(con[,1]=="Matrix for phylo"),2])
}else{
	textEntryKIB <- tclVar("")
}
textEntryWidgetKIB <- tkentry(epi,width=10,textvariable=textEntryKIB)
tkconfigure(textEntryWidgetKIB, state="disabled")
tkgrid(pElabel1KIB, textEntryWidgetKIB)

pElabelTextH42B <- tclVar("Window size: ")
pElabelH42B <- tklabel(epi,text=tclvalue(pElabelTextH42B))
if (length(con[which(con[,1]=="window_size"),2]) > 0){
	textEntryH42B <- tclVar(con[which(con[,1]=="window_size"),2])
}else{
	textEntryH42B <- tclVar("9")
}
textEntryWidgetH42B <- tkentry(epi,width=10,textvariable=textEntryH42B)
tkgrid(pElabelH42B, textEntryWidgetH42B)


OK.bute1 <- tkbutton(epi,text="Start",command=function()calculate_pos_epi(tbn, cbValue01, cbValue_0e, textEntryC, textEntryD, textEntryE, textEntryF, textEntryG, textEntryH, textEntryH1, textEntryH2, textEntryH12, textEntryHB, textEntryH1B, textEntryH2B, textEntryH12B, textEntryIB, textEntryKIB, txt, as.character(tclvalue(icbValue01)), textEntryH42B ))
tkgrid(OK.bute1)

tkgrid.configure(button.widget_epi1, column=7)
tkgrid.configure(button.widget_epi2, column=7)
tkgrid.configure(button.widget_epi3, column=7)
tkgrid.configure(button.widget_epi4, column=7)
tkgrid.configure(pElabel2, columnspan=5)
tkgrid.configure(pElabel2A, columnspan=5)
tkgrid.configure(pElabel2B, columnspan=5)
tkgrid.configure(textEntryWidgetC, columnspan=5)
tkgrid.configure(textEntryWidgetD, columnspan=5)
tkgrid.configure(textEntryWidgetE, columnspan=5)
tkgrid.configure(textEntryWidgetF, columnspan=5)
tkgrid.configure(textEntryWidgetG, columnspan=5)
tkgrid.configure(ddb, columnspan=5)

#-----------------------------------------------------------
lb1 <- tklabel(epist, text=paste("Create two-side plot"), background="#FFCC99")
tkgrid (lb1, sticky="snew", columnspan=8)

Bbutton.widget_epi <- tkbutton(epist,text="Select csv file",command=function()getcsvfile(BpElabelText1, 7, 1))
BpElabelText <- tclVar("Selected epitope results: ")
BpElabel1 <- tklabel(epist,text=tclvalue(BpElabelText))
tkconfigure(BpElabel1,textvariable=BpElabelText)

BpElabelText1 <- tclVar("")
BpElabel2 <- tklabel(epist,text=tclvalue(BpElabelText1))
tkconfigure(BpElabel2,textvariable=BpElabelText1)
tkgrid(BpElabel1, BpElabel2, Bbutton.widget_epi)

BpElabelTextA <- tclVar("Csv seperator: ")
BpElabel1A <- tklabel(epist,text=tclvalue(BpElabelTextA))
if (length(con[which(con[,1]=="csv seperator"),2]) > 0){
	BtextEntryA <- tclVar(con[which(con[,1]=="csv seperator"),2])
}else{
	BtextEntryA <- tclVar("")
}
BtextEntryWidgetA <- tkentry(epist,width=5,textvariable=BtextEntryA)
tkgrid(BpElabel1A, BtextEntryWidgetA)

BpElabelTextB <- tclVar("Number of cases: ")
BpElabel1B <- tklabel(epist,text=tclvalue(BpElabelTextB))
if (length(con[which(con[,1]=="number of cases"),2]) > 0){
	BtextEntryB <- tclVar(con[which(con[,1]=="number of cases"),2])
}else{
	BtextEntryB <- tclVar("")
}
BtextEntryWidgetB <- tkentry(epist,width=5,textvariable=BtextEntryB)
tkgrid(BpElabel1B, BtextEntryWidgetB)

BpElabelTextC <- tclVar("First column of odds-ratio: ")
BpElabel1C <- tklabel(epist,text=tclvalue(BpElabelTextC))
if (length(con[which(con[,1]=="position of odds.ratio"),2]) > 0){
	BtextEntryC <- tclVar(con[which(con[,1]=="position of odds.ratio"),2])
}else{
	BtextEntryC <- tclVar("")
}
BtextEntryWidgetC <- tkentry(epist,width=5,textvariable=BtextEntryC)
tkgrid(BpElabel1C, BtextEntryWidgetC)

BpElabelTextD <- tclVar("First column of p-values: ")
BpElabel1D <- tklabel(epist,text=tclvalue(BpElabelTextD))
if (length(con[which(con[,1]=="position of p.values"),2]) > 0){
	BtextEntryD <- tclVar(con[which(con[,1]=="position of p.values"),2])
}else{
	BtextEntryD <- tclVar("")
}
BtextEntryWidgetD <- tkentry(epist,width=5,textvariable=BtextEntryD)
tkgrid(BpElabel1D, BtextEntryWidgetD)

BpElabelTextE <- tclVar("Column for name for y-axis label: ")
BpElabel1E <- tklabel(epist,text=tclvalue(BpElabelTextE))
if (length(con[which(con[,1]=="name for the y-axis label"),2]) > 0){
	BtextEntryE <- tclVar(con[which(con[,1]=="name for the y-axis label"),2])
}else{
	BtextEntryE <- tclVar("")
}
BtextEntryWidgetE <- tkentry(epist,width=5,textvariable=BtextEntryE)
tkgrid(BpElabel1E, BtextEntryWidgetE)

BpElabelTextF <- tclVar("Number of columns for one case: ")
BpElabel1F <- tklabel(epist,text=tclvalue(BpElabelTextF))
if (length(con[which(con[,1]=="frequency"),2]) > 0){
	BtextEntryF <- tclVar(con[which(con[,1]=="frequency"),2])
}else{
	BtextEntryF <- tclVar("")
}
BtextEntryWidgetF <- tkentry(epist,width=5,textvariable=BtextEntryF)
tkgrid(BpElabel1F, BtextEntryWidgetF)

BpElabelTextG <- tclVar("First column of list of aa: ")
BpElabel1G <- tklabel(epist,text=tclvalue(BpElabelTextG))
if (length(con[which(con[,1]=="position of amino acid"),2]) > 0){
	BtextEntryG <- tclVar(con[which(con[,1]=="position of amino acid"),2])
}else{
	BtextEntryG <- tclVar("")
}
BtextEntryWidgetG <- tkentry(epist,width=5,textvariable=BtextEntryG)
tkgrid(BpElabel1G, BtextEntryWidgetG)

BpElabelTextH <- tclVar("Maximum value of y-axis: ")
BpElabel1H <- tklabel(epist,text=tclvalue(BpElabelTextH))
if (length(con[which(con[,1]=="maximum value of y-axis"),2]) > 0){
	BtextEntryH <- tclVar(con[which(con[,1]=="maximum value of y-axis"),2])
}else{
	BtextEntryH <- tclVar("")
}
BtextEntryWidgetH <- tkentry(epist,width=5,textvariable=BtextEntryH)
tkgrid(BpElabel1H, BtextEntryWidgetH)

BpElabelTextI <- tclVar("Ticks of y-axis: ")
BpElabel1I <- tklabel(epist,text=tclvalue(BpElabelTextI))
if (length(con[which(con[,1]=="intervall of y-axis"),2]) > 0){
	BtextEntryI <- tclVar(con[which(con[,1]=="intervall of y-axis"),2])
}else{
	BtextEntryI <- tclVar("")
}
BtextEntryWidgetI <- tkentry(epist,width=5,textvariable=BtextEntryI)
tkgrid(BpElabel1I, BtextEntryWidgetI)

BpElabelTextJ <- tclVar("First column for additional color information: ")
BpElabel1J <- tklabel(epist,text=tclvalue(BpElabelTextJ))
if (length(con[which(con[,1]=="add color information"),2]) > 0){
	BtextEntryJ <- tclVar(con[which(con[,1]=="add color information"),2])
}else{
	BtextEntryJ <- tclVar("")
}
BtextEntryWidgetJ <- tkentry(epist,width=5,textvariable=BtextEntryJ)
tkgrid(BpElabel1J, BtextEntryWidgetJ)

BOK.bute2 <- tkbutton(epist,text="Start",command=function()create_sequence_graphic_G(tbn, BtextEntryA, BtextEntryB, BtextEntryC, BtextEntryD, BtextEntryE, BtextEntryF, BtextEntryG, BtextEntryH, BtextEntryI, txt, BtextEntryJ))
tkgrid(BOK.bute2)


#--------------------------------------------------------------------------------------------------------------------------------
tkgrid.rowconfigure(otbn2,0,weight=1, minsize=100)
tkgrid.columnconfigure(otbn2,0,weight=10, minsize=200)
tkgrid.columnconfigure(otbn2,1,weight=3, minsize=100)

tbn2 <- ttkframe(otbn2)
helptbn2 <- ttkframe(otbn2)
tkgrid(tbn2, helptbn2, sticky="snew")
tkgrid.columnconfigure(tbn2,0,weight=10)
tkgrid.columnconfigure(helptbn2,0,weight=1)

insidetbn2 <- tkframe(tbn2)
tkpack(insidetbn2,fill='both',expand='yes')

scr2 <- tkscrollbar(insidetbn2, command=function(...)tkyview(txt2,...),bg='white')
txt2 <- tktext(insidetbn2,height=1, width=99,bg='grey')
xscr2 <- tkscrollbar(insidetbn2, command=function(...)tkxview(txt2,...),orient='horiz')
tkconfigure(txt2, yscrollcommand=function(...)tkset(scr2,...), xscrollcommand=function(...)tkset(xscr2,...), state="disabled")
tkpack(scr2,side='right',fill='y')
tkpack(txt2,expand='yes',fill='both')
tkpack(xscr2, side="bottom", fill="x")
sf <- tkcanvas(txt2, background="white")
tkwindow.create(txt2, 'end', window=sf)

cb1 <- tkcheckbutton(sf,command=function()check_allels(cbValu <- as.character(tclvalue(cbValuee)), cMtextEntryWidgetC, button.widget_co2, textEntryWidgetHC, textEntryWidgetH1C, textEntryWidgetH2C, textEntryWidgetH12C, textEntryWidgetHBC, textEntryWidgetH1BC, textEntryWidgetH2BC, textEntryWidgetH12BC), background="#FF6600")
if (length(con[which(con[,1]=="With Allels"),2]) > 0){
	cbValuee <- tclVar(con[which(con[,1]=="With Allels"),2])
}else{
	cbValuee <- tclVar("0")
}
tkconfigure(cb1,variable=cbValuee)
tkgrid(tklabel(sf,text="Analysis with more than one feature in data (aka HLA allels):", background="#FF6600"),cb1, sticky="snew", columnspan=1)

calc <- tkframe(sf,relief="groove",borderwidth=2)
sm <- tkframe(sf,relief="groove",borderwidth=2)

tkgrid(calc, sticky="snew")
tkgrid(sm, sticky="snew")

heading1 <- tklabel(calc,text="Discovering associations between mutation pairs and features", background="#FFCC99")
tkgrid(heading1, sticky="snew", columnspan=8)

button.widget_co1 <- tkbutton(calc,text="Select FASTA file",command=function()getFASTAfile(cMlabelText1, 2, 1))

cMlabelText <- tclVar("Selected sequences: ")
cMlabel1 <- tklabel(calc,text=tclvalue(cMlabelText))
tkconfigure(cMlabel1,textvariable=cMlabelText)

cMlabelText1 <- tclVar("")
cMlabel2 <- tklabel(calc,text=tclvalue(cMlabelText1))
tkconfigure(cMlabel2,textvariable=cMlabelText1)
tkgrid(cMlabel1, cMlabel2, button.widget_co1)

button.widget_co2 <- tkbutton(calc,text="Select FASTA file",command=function()getFASTAfile(cMlabelText1A, 2, 2))
tkconfigure(button.widget_co2, state="normal")

cMlabelTextA <- tclVar("Selected consensus: ")
cMlabel1A <- tklabel(calc,text=tclvalue(cMlabelTextA))
tkconfigure(cMlabel1A,textvariable=cMlabelTextA)

cMlabelText1A <- tclVar("")
cMlabel2A <- tklabel(calc,text=tclvalue(cMlabelText1A))
tkconfigure(cMlabel2A,textvariable=cMlabelText1A)
tkgrid(cMlabel1A, cMlabel2A, button.widget_co2)

cMlabelTextC <- tclVar("Minimal number of members: ")
cMlabel1C <- tklabel(calc,text=tclvalue(cMlabelTextC))
if (length(con[which(con[,1]=="Number of patients threshold_cm"),2]) > 0){
	cMtextEntryC <- tclVar(con[which(con[,1]=="Number of patients threshold_cm"),2])
}else{
	cMtextEntryC <- tclVar("")
}
cMtextEntryWidgetC <- tkentry(calc,width=20,textvariable=cMtextEntryC)
tkconfigure(cMtextEntryWidgetC, state="disabled")
tkgrid(cMlabel1C, cMtextEntryWidgetC)

cMlabelTextD <- tclVar("Level for significance: ")
cMlabel1D <- tklabel(calc,text=tclvalue(cMlabelTextD))
if (length(con[which(con[,1]=="level for significance"),2]) > 0){
	cMtextEntryD <- tclVar(con[which(con[,1]=="level for significance"),2])
}else{
	cMtextEntryD <- tclVar("")
}
cMtextEntryWidgetD <- tkentry(calc,width=20,textvariable=cMtextEntryD)
tkgrid(cMlabel1D, cMtextEntryWidgetD)

cbm <- tkcheckbutton(calc)
if (length(con[which(con[,1]=="More than one core"),2]) > 0){
	cbValuem <- tclVar(con[which(con[,1]=="More than one core"),2])
}else{
	cbValuem <- tclVar("0")
}
tkconfigure(cbm,variable=cbValuem)
tkgrid(tklabel(calc,text="Use more than one core?"),cbm)

cMlabelTextHC <- tclVar("Position of feature A: ")
cMlabel1HC <- tklabel(calc,text=tclvalue(cMlabelTextHC))
if (length(con[which(con[,1]=="cmHLA-A1"),2]) > 0){
	textEntryHC <- tclVar(con[which(con[,1]=="cmHLA-A1"),2])
}else{
	textEntryHC <- tclVar("")
}
textEntryWidgetHC <- tkentry(calc,width=3,textvariable=textEntryHC)
cMlabelTextH1C <- tclVar("-")
cMlabel1H1C <- tklabel(calc,text=tclvalue(cMlabelTextH1C))
if (length(con[which(con[,1]=="cmHLA-A2"),2]) > 0){
	textEntryH1C <- tclVar(con[which(con[,1]=="cmHLA-A2"),2])
}else{
	textEntryH1C <- tclVar("")
}
textEntryWidgetH1C <- tkentry(calc,width=3,textvariable=textEntryH1C)
if (length(con[which(con[,1]=="cmHLA-A3"),2]) > 0){
	textEntryH2C <- tclVar(con[which(con[,1]=="cmHLA-A3"),2])
}else{
	textEntryH2C <- tclVar("")
}
textEntryWidgetH2C <- tkentry(calc,width=3,textvariable=textEntryH2C)
cMlabelTextH12C <- tclVar("-")
cMlabel1H12C <- tklabel(calc,text=tclvalue(cMlabelTextH12C))
if (length(con[which(con[,1]=="cmHLA-A4"),2]) > 0){
	textEntryH12C <- tclVar(con[which(con[,1]=="cmHLA-A4"),2])
}else{
	textEntryH12C <- tclVar("")
}
textEntryWidgetH12C <- tkentry(calc,width=3,textvariable=textEntryH12C)
tkgrid(cMlabel1HC, textEntryWidgetHC, cMlabel1H1C, textEntryWidgetH1C, textEntryWidgetH2C, cMlabel1H12C, textEntryWidgetH12C)
tkconfigure(textEntryWidgetHC, state="disabled")
tkconfigure(textEntryWidgetH1C, state="disabled")
tkconfigure(textEntryWidgetH2C, state="disabled")
tkconfigure(textEntryWidgetH12C, state="disabled")

cMlabelTextHBC <- tclVar("Position of feature B: ")
cMlabel1HBC <- tklabel(calc,text=tclvalue(cMlabelTextHBC))
if (length(con[which(con[,1]=="cmHLA-B1"),2]) > 0){
	textEntryHBC <- tclVar(con[which(con[,1]=="cmHLA-B1"),2])
}else{
	textEntryHBC <- tclVar("")
}
textEntryWidgetHBC <- tkentry(calc,width=3,textvariable=textEntryHBC)
cMlabelTextH1BC <- tclVar("-")
cMlabel1H1BC <- tklabel(calc,text=tclvalue(cMlabelTextH1BC))
if (length(con[which(con[,1]=="cmHLA-B2"),2]) > 0){
	textEntryH1BC <- tclVar(con[which(con[,1]=="cmHLA-B2"),2])
}else{
	textEntryH1BC <- tclVar("")
}
textEntryWidgetH1BC <- tkentry(calc,width=3,textvariable=textEntryH1BC)
if (length(con[which(con[,1]=="cmHLA-B3"),2]) > 0){
	textEntryH2BC <- tclVar(con[which(con[,1]=="cmHLA-B3"),2])
}else{
	textEntryH2BC <- tclVar("")
}
textEntryWidgetH2BC <- tkentry(calc,width=3,textvariable=textEntryH2BC)
cMlabelTextH12BC <- tclVar("-")
cMlabel1H12BC <- tklabel(calc,text=tclvalue(cMlabelTextH12BC))
if (length(con[which(con[,1]=="cmHLA-B4"),2]) > 0){
	textEntryH12BC <- tclVar(con[which(con[,1]=="cmHLA-B4"),2])
}else{
	textEntryH12BC <- tclVar("")
}
textEntryWidgetH12BC <- tkentry(calc,width=3,textvariable=textEntryH12BC)
tkgrid(cMlabel1HBC, textEntryWidgetHBC, cMlabel1H1BC, textEntryWidgetH1BC, textEntryWidgetH2BC, cMlabel1H12BC, textEntryWidgetH12BC)
tkconfigure(textEntryWidgetHBC, state="disabled")
tkconfigure(textEntryWidgetH1BC, state="disabled")
tkconfigure(textEntryWidgetH2BC, state="disabled")
tkconfigure(textEntryWidgetH12BC, state="disabled")

cmlabelTextIB <- tclVar("P-value correction: ")
cmlabel1IB <- tklabel(calc,text=tclvalue(cmlabelTextIB))
ddbc <- tkframe(calc)
cmtextEntryIB <- tclVar()
dropdownList(ddbc, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), cmtextEntryIB, 15, "none")
tkgrid(cmlabel1IB, ddbc)

OK.but21 <- tkbutton(calc,text="Start",command=function()calculate_co_mut(cbVal <- as.character(tclvalue(cbValuee)), tbn2, txt, cMtextEntryC, cMtextEntryD, cmtextEntryIB, textEntryHC, textEntryH1C, textEntryH2C, textEntryH12C, textEntryHBC, textEntryH1BC, textEntryH2BC, textEntryH12BC, cmtextEntryIB, cbValue01))
tkgrid(OK.but21)

heading4 <- tklabel(sm,text="Discovering associations between mutation tuples and features", background="#FFCC99")
tkgrid(heading4, sticky="snew", columnspan=8)

button.widget_sm1 <- tkbutton(sm,text="Select original FASTA file",command=function()getFASTAfile(smlabelText1, 8, 1))
smlabelText <- tclVar("Original FASTA file: ")
smlabel1 <- tklabel(sm,text=tclvalue(smlabelText))
tkconfigure(smlabel1,textvariable=smlabelText)
smlabelText1 <- tclVar("")
smlabel12 <- tklabel(sm,text=tclvalue(smlabelText1))
tkconfigure(smlabel12,textvariable=smlabelText1)
tkgrid(smlabel1, smlabel12, button.widget_sm1)


button.widget_sm2 <- tkbutton(sm,text="Select csv result file",command=function()getcsvfile(smlabelText3, 8, 2))
smlabelText2 <- tclVar("Results from \"point mutations\": ")
smlabel2 <- tklabel(sm,text=tclvalue(smlabelText2))
tkconfigure(smlabel2,textvariable=smlabelText2)
smlabelText3 <- tclVar("")
smlabel3 <- tklabel(sm,text=tclvalue(smlabelText3))
tkconfigure(smlabel3,textvariable=smlabelText3)
tkgrid(smlabel2, smlabel3, button.widget_sm2)

smlabelTextA <- tclVar("Threshold for position selection: ")
smlabel1A <- tklabel(sm,text=tclvalue(smlabelTextA))
if (length(con[which(con[,1]=="threshold"),2]) > 0){
	smtextEntryA <- tclVar(con[which(con[,1]=="threshold"),2])
}else{
	smtextEntryA <- tclVar("")
}
smtextEntryWidgetA <- tkentry(sm,width=20,textvariable=smtextEntryA)
tkgrid(smlabel1A, smtextEntryWidgetA)

smlabelTextB <- tclVar("Minimum number of elements in tuple: ")
smlabel1B <- tklabel(sm,text=tclvalue(smlabelTextB))
if (length(con[which(con[,1]=="min_number_of_ele_in_tupel"),2]) > 0){
	smtextEntryB <- tclVar(con[which(con[,1]=="min_number_of_ele_in_tupel"),2])
}else{
	smtextEntryB <- tclVar("")
}
smtextEntryWidgetB <- tkentry(sm,width=20,textvariable=smtextEntryB)
tkgrid(smlabel1B, smtextEntryWidgetB)

smlabelTextC <- tclVar("Maximum number of elements in tuple: ")
smlabel1C <- tklabel(sm,text=tclvalue(smlabelTextC))
if (length(con[which(con[,1]=="max_number_of_ele_in_tuple"),2]) > 0){
	smtextEntryC <- tclVar(con[which(con[,1]=="max_number_of_ele_in_tuple"),2])
}else{
	smtextEntryC <- tclVar("")
}
smtextEntryWidgetC <- tkentry(sm,width=20,textvariable=smtextEntryC)
tkgrid(smlabel1C, smtextEntryWidgetC)

smlabelTextD <- tclVar("Column in which the identifier is noted: ")
smlabel1D <- tklabel(sm,text=tclvalue(smlabelTextD))
if (length(con[which(con[,1]=="column"),2]) > 0){
	smtextEntryD <- tclVar(con[which(con[,1]=="column"),2])
}else{
	smtextEntryD <- tclVar("")
}
smtextEntryWidgetD <- tkentry(sm,width=20,textvariable=smtextEntryD)
tkgrid(smlabel1D, smtextEntryWidgetD)

smlabelTextE <- tclVar("Column of p-values for this identifier: ")
smlabel1E <- tklabel(sm,text=tclvalue(smlabelTextE))
if (length(con[which(con[,1]=="column of values"),2]) > 0){
	smtextEntryE <- tclVar(con[which(con[,1]=="column of values"),2])
}else{
	smtextEntryE <- tclVar("")
}
smtextEntryWidgetE <- tkentry(sm,width=20,textvariable=smtextEntryE)
tkgrid(smlabel1E, smtextEntryWidgetE)

smlabelTextF <- tclVar("Column of aa/nt for this identifier: ")
smlabel1F <- tklabel(sm,text=tclvalue(smlabelTextF))
if (length(con[which(con[,1]=="column of aas"),2]) > 0){
	smtextEntryF <- tclVar(con[which(con[,1]=="column of aas"),2])
}else{
	smtextEntryF <- tclVar("")
}
smtextEntryWidgetF <- tkentry(sm,width=20,textvariable=smtextEntryF)
tkgrid(smlabel1F, smtextEntryWidgetF)

smlabelTextH <- tclVar("Position of feature A: ")
smlabel1H <- tklabel(sm,text=tclvalue(smlabelTextH))
if (length(con[which(con[,1]=="smHLA-A1"),2]) > 0){
	smtextEntryH <- tclVar(con[which(con[,1]=="smHLA-A1"),2])
}else{
	smtextEntryH <- tclVar("")
}
smtextEntryWidgetH <- tkentry(sm,width=3,textvariable=smtextEntryH)
smlabelTextH1 <- tclVar("-")
smlabel1H1 <- tklabel(sm,text=tclvalue(smlabelTextH1))
if (length(con[which(con[,1]=="smHLA-A2"),2]) > 0){
	smtextEntryH1 <- tclVar(con[which(con[,1]=="smHLA-A2"),2])
}else{
	smtextEntryH1 <- tclVar("")
}
smtextEntryWidgetH1 <- tkentry(sm,width=3,textvariable=smtextEntryH1)
if (length(con[which(con[,1]=="smHLA-A3"),2]) > 0){
	smtextEntryH2 <- tclVar(con[which(con[,1]=="smHLA-A3"),2])
}else{
	smtextEntryH2 <- tclVar("")
}
smtextEntryWidgetH2 <- tkentry(sm,width=3,textvariable=smtextEntryH2)
smlabelTextH12 <- tclVar("-")
smlabel1H12 <- tklabel(sm,text=tclvalue(smlabelTextH12))
if (length(con[which(con[,1]=="smHLA-A4"),2]) > 0){
	smtextEntryH12 <- tclVar(con[which(con[,1]=="smHLA-A4"),2])
}else{
	smtextEntryH12 <- tclVar("")
}
smtextEntryWidgetH12 <- tkentry(sm,width=3,textvariable=textEntryH12)
tkgrid(smlabel1H, smtextEntryWidgetH, smlabel1H1, smtextEntryWidgetH1, smtextEntryWidgetH2, smlabel1H12, smtextEntryWidgetH12)

smlabelTextHB <- tclVar("Position of feature B: ")
smlabel1HB <- tklabel(sm,text=tclvalue(smlabelTextHB))
if (length(con[which(con[,1]=="smHLA-B1"),2]) > 0){
	smtextEntryHB <- tclVar(con[which(con[,1]=="smHLA-B1"),2])
}else{
	smtextEntryHB <- tclVar("")
}
smtextEntryWidgetHB <- tkentry(sm,width=3,textvariable=smtextEntryHB)
smlabelTextH1B <- tclVar("-")
smlabel1H1B <- tklabel(sm,text=tclvalue(smlabelTextH1B))
if (length(con[which(con[,1]=="smHLA-B2"),2]) > 0){
	smtextEntryH1B <- tclVar(con[which(con[,1]=="smHLA-B2"),2])
}else{
	smtextEntryH1B <- tclVar("")
}
smtextEntryWidgetH1B <- tkentry(sm,width=3,textvariable=smtextEntryH1B)
if (length(con[which(con[,1]=="smHLA-B3"),2]) > 0){
	smtextEntryH2B <- tclVar(con[which(con[,1]=="smHLA-B3"),2])
}else{
	smtextEntryH2B <- tclVar("")
}
smtextEntryWidgetH2B <- tkentry(sm,width=3,textvariable=smtextEntryH2B)
smlabelTextH12B <- tclVar("-")
smlabel1H12B <- tklabel(sm,text=tclvalue(smlabelTextH12B))
if (length(con[which(con[,1]=="smHLA-B4"),2]) > 0){
	smtextEntryH12B <- tclVar(con[which(con[,1]=="smHLA-B4"),2])
}else{
	smtextEntryH12B <- tclVar("")
}
smtextEntryWidgetH12B <- tkentry(sm,width=3,textvariable=smtextEntryH12B)
tkgrid(smlabel1HB, smtextEntryWidgetHB, smlabel1H1B, smtextEntryWidgetH1B, smtextEntryWidgetH2B, smlabel1H12B, smtextEntryWidgetH12B)

smlabelTextX <- tclVar("Identifier (if only one feature): ")
smlabelX <- tklabel(sm,text=tclvalue(smlabelTextX))
if (length(con[which(con[,1]=="Identifier"),2]) > 0){
	smtextEntryidet <- tclVar(con[which(con[,1]=="Identifier"),2])
}else{
	smtextEntryidet <- tclVar("")
}
smtextEntryWidgetX <- tkentry(sm,width=3,textvariable=smtextEntryidet)
tkgrid(smlabelX, smtextEntryWidgetX)
tkconfigure(smtextEntryWidgetX, state="disabled")

OK.but4 <- tkbutton(sm,text="Start",command=function()shared_mutations(tbn2, txt, smtextEntryA, smtextEntryB, smtextEntryC, smtextEntryD, smtextEntryE, smtextEntryF, smtextEntryH, smtextEntryH1, smtextEntryH2, smtextEntryH12, smtextEntryHB, smtextEntryH1B, smtextEntryH2B, smtextEntryH12B, as.character(tclvalue(icbValue01)), smtextEntryidet ))
tkgrid(OK.but4)

tkgrid.configure(button.widget_sm1, column=7)
tkgrid.configure(button.widget_sm2, column=7)
tkgrid.configure(smlabel12, columnspan=5)
tkgrid.configure(smlabel3, columnspan=5)
tkgrid.configure(smtextEntryWidgetA, columnspan=5)
tkgrid.configure(smtextEntryWidgetB, columnspan=5)
tkgrid.configure(smtextEntryWidgetC, columnspan=5)
tkgrid.configure(smtextEntryWidgetD, columnspan=5)
tkgrid.configure(smtextEntryWidgetE, columnspan=5)
tkgrid.configure(smtextEntryWidgetF, columnspan=5)
#--------------------------------------------------------------------------------------------------------------------------------
tkgrid.rowconfigure(otbn5,0,weight=1, minsize=100)
tkgrid.columnconfigure(otbn5,0,weight=10, minsize=200)
tkgrid.columnconfigure(otbn5,1,weight=3, minsize=100)

tbn5 <- ttkframe(otbn5)
helptbn5 <- ttkframe(otbn5)
tkgrid(tbn5, helptbn5, sticky="snew")
tkgrid.columnconfigure(tbn5,0,weight=10)
tkgrid.columnconfigure(helptbn5,0,weight=1)

scr3 <- tkscrollbar(tbn5, command=function(...)tkyview(txt3,...),bg='white')
txt3 <- tktext(tbn5,height=1,width=55,bg='grey') 
xscr3 <- tkscrollbar(tbn5, command=function(...)tkxview(txt3,...),orient='horiz')
tkconfigure(txt3, yscrollcommand=function(...)tkset(scr3,...), xscrollcommand=function(...)tkset(xscr3,...), state="disabled")
tkpack(scr3,side='right',fill='y')
tkpack(txt3,expand='yes',fill='both')
tkpack(xscr3, side="bottom", fill="x")
sf3 <- tkcanvas(txt3, background="white")
tkwindow.create(txt3, 'end', window=sf3)

cb <- tkcheckbutton(sf3,command=function()check_allels(cbValu <- as.character(tclvalue(cbValue)), cMtextEntryWidgetC, button.widget_co2, textEntryWidgetHC, textEntryWidgetH1C, textEntryWidgetH2C, textEntryWidgetH12C, textEntryWidgetHBC, textEntryWidgetH1BC, textEntryWidgetH2BC, textEntryWidgetH12BC), background="#FF6600")
if (length(con[which(con[,1]=="cbWith Allels"),2]) > 0){
	cbValue <- tclVar(con[which(con[,1]=="cbWith Allels"),2])
}else{
	cbValue <- tclVar("0")
}
tkconfigure(cb,variable=cbValue)
tkgrid(tklabel(sf3,text="Analysis with more than one feature in data (aka HLA allels):", background="#FF6600"),cb, sticky="snew", columnspan=1)

graphicals <- tkframe(sf3,relief="groove",borderwidth=2)
graphicals2 <- tkframe(sf3,relief="groove",borderwidth=2)

tkgrid(graphicals, sticky="snew")
tkgrid(graphicals2, sticky="snew")
#Co mutation graphics

heading2 <- tklabel(graphicals,text="Create graphics for tuples, one page per identifier", background="#FFCC99")
tkgrid(heading2, sticky="snew", columnspan=3)

button.widget_com <- tkbutton(graphicals,text="Select csv file",command=function()getcsvfile(cMlabelText1B, 5, 1))

cMlabelTextB <- tclVar("Selected result from tuple analysis: ")
cMlabel1B <- tklabel(graphicals,text=tclvalue(cMlabelTextB))
tkconfigure(cMlabel1B,textvariable=cMlabelTextB)

cMlabelText1B <- tclVar("")
cMlabel2B <- tklabel(graphicals,text=tclvalue(cMlabelText1B))
tkconfigure(cMlabel2B,textvariable=cMlabelText1B)
tkgrid(cMlabel1B, cMlabel2B, button.widget_com)

cMlabelTextE <- tclVar("Level for significance: ")
cMlabel1E <- tklabel(graphicals,text=tclvalue(cMlabelTextE))
if (length(con[which(con[,1]=="cMlevel for significance"),2]) > 0){
	cMtextEntryE <- tclVar(con[which(con[,1]=="cMlevel for significance"),2])
}else{
	cMtextEntryE <- tclVar("")
}
cMtextEntryWidgetE <- tkentry(graphicals,width=20,textvariable=cMtextEntryE)
tkgrid(cMlabel1E, cMtextEntryWidgetE)

OK.but22 <- tkbutton(graphicals,text="Start",command=function()calculate_co_mut_graphics(cbVal2 <- as.character(tclvalue(cbValue)), tbn2, cMtextEntryE))
tkgrid(OK.but22)

#double co mutation graphics

heading4 <- tklabel(graphicals2,text="Create tartan plot", background="#FFCC99")
tkgrid(heading4, sticky="snew", columnspan=3)

Tbutton.widget_com <- tkbutton(graphicals2,text="Select csv file",command=function()getcsvfile(TcMlabelText1B, 9, 1))
Tbutton.widget_com2 <- tkbutton(graphicals2,text="Select csv file",command=function()getcsvfile(TcMlabelText1B2, 9, 2))
Tbutton.widget_com3 <- tkbutton(graphicals2,text="Select csv file",command=function()getcsvfile(TcMlabelText1B3, 9, 3))

TcMlabelTextB <- tclVar("Selected 1st result from co-mutation: ")
TcMlabel1B <- tklabel(graphicals2,text=tclvalue(TcMlabelTextB))
tkconfigure(TcMlabel1B,textvariable=TcMlabelTextB)

TcMlabelText1B <- tclVar("")
TcMlabel2B <- tklabel(graphicals2,text=tclvalue(TcMlabelText1B))
tkconfigure(TcMlabel2B,textvariable=TcMlabelText1B)
tkgrid(TcMlabel1B, TcMlabel2B, Tbutton.widget_com)

TcMlabelTextB2 <- tclVar("Selected 2nd result from co-mutation: ")
TcMlabel1B2 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextB2))
tkconfigure(TcMlabel1B2,textvariable=TcMlabelTextB2)

TcMlabelText1B2 <- tclVar("")
TcMlabel2B2 <- tklabel(graphicals2,text=tclvalue(TcMlabelText1B2))
tkconfigure(TcMlabel2B2,textvariable=TcMlabelText1B2)
tkgrid(TcMlabel1B2, TcMlabel2B2, Tbutton.widget_com2)

TcMlabelTextB3 <- tclVar("Optional selected distance matrix: ")
TcMlabel1B3 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextB3))
tkconfigure(TcMlabel1B3,textvariable=TcMlabelTextB3)

TcMlabelText1B3 <- tclVar("")
TcMlabel2B3 <- tklabel(graphicals2,text=tclvalue(TcMlabelText1B3))
tkconfigure(TcMlabel2B3,textvariable=TcMlabelText1B3)
tkgrid(TcMlabel1B3, TcMlabel2B3, Tbutton.widget_com3)

TcMlabelTextE <- tclVar("Optical space between blocks: ")
TcMlabel1E <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE))
if (length(con[which(con[,1]=="space between blocks"),2]) > 0){
	TcMtextEntryE <- tclVar(con[which(con[,1]=="space between blocks"),2])
}else{
	TcMtextEntryE <- tclVar("")
}
TcMtextEntryWidgetE <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE)
tkgrid(TcMlabel1E, TcMtextEntryWidgetE)

TcMlabelTextE1 <- tclVar("Colors of the plot: ")
TcMlabel1E1 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE1))
if (length(con[which(con[,1]=="colors of the plot"),2]) > 0){
	TcMtextEntryE1 <- tclVar(con[which(con[,1]=="colors of the plot"),2])
}else{
	TcMtextEntryE1 <- tclVar("")
}
TcMtextEntryWidgetE1 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE1)
tkgrid(TcMlabel1E1, TcMtextEntryWidgetE1)

TcMlabelTextE2 <- tclVar("x/y positions where labels should be added: ")
TcMlabel1E2 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE2))
if (length(con[which(con[,1]=="name positions"),2]) > 0){
	TcMtextEntryE2 <- tclVar(con[which(con[,1]=="name positions"),2])
}else{
	TcMtextEntryE2 <- tclVar("")
}
TcMtextEntryWidgetE2 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE2)
tkgrid(TcMlabel1E2, TcMtextEntryWidgetE2)

TcMlabelTextE3 <- tclVar("Labels to add in above positions: ")
TcMlabel1E3 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE3))
if (length(con[which(con[,1]=="names"),2]) > 0){
	TcMtextEntryE3 <- tclVar(con[which(con[,1]=="names"),2])
}else{
	TcMtextEntryE3 <- tclVar("")
}
TcMtextEntryWidgetE3 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE3)
tkgrid(TcMlabel1E3, TcMtextEntryWidgetE3)

TcMlabelTextE4 <- tclVar("Ticks: ")
TcMlabel1E4 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE4))
if (length(con[which(con[,1]=="ticks"),2]) > 0){
	TcMtextEntryE4 <- tclVar(con[which(con[,1]=="ticks"),2])
}else{
	TcMtextEntryE4 <- tclVar("")
}
TcMtextEntryWidgetE4 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE4)
tkgrid(TcMlabel1E4, TcMtextEntryWidgetE4)

TcMlabelTextE5 <- tclVar("Column of 1st sequence position in 1st file: ")
TcMlabel1E5 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE5))
if (length(con[which(con[,1]=="column of first position in first file"),2]) > 0){
	TcMtextEntryE5 <- tclVar(con[which(con[,1]=="column of first position in first file"),2])
}else{
	TcMtextEntryE5 <- tclVar("")
}
TcMtextEntryWidgetE5 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE5)
tkgrid(TcMlabel1E5, TcMtextEntryWidgetE5)

TcMlabelTextE6 <- tclVar("Column of 2nd sequence position in 1st file: ")
TcMlabel1E6 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE6))
if (length(con[which(con[,1]=="column of second position in first file"),2]) > 0){
	TcMtextEntryE6 <- tclVar(con[which(con[,1]=="column of second position in first file"),2])
}else{
	TcMtextEntryE6 <- tclVar("")
}
TcMtextEntryWidgetE6 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE6)
tkgrid(TcMlabel1E6, TcMtextEntryWidgetE6)

TcMlabelTextE7 <- tclVar("Column of p-values in 1st file: ")
TcMlabel1E7 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE7))
if (length(con[which(con[,1]=="column of values in first file"),2]) > 0){
	TcMtextEntryE7 <- tclVar(con[which(con[,1]=="column of values in first file"),2])
}else{
	TcMtextEntryE7 <- tclVar("")
}
TcMtextEntryWidgetE7 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE7)
tkgrid(TcMlabel1E7, TcMtextEntryWidgetE7)

TcMlabelTextE8 <- tclVar("Column of 1st sequence position in 2nd file: ")
TcMlabel1E8 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE8))
if (length(con[which(con[,1]=="column of first position in second file"),2]) > 0){
	TcMtextEntryE8 <- tclVar(con[which(con[,1]=="column of first position in second file"),2])
}else{
	TcMtextEntryE8 <- tclVar("")
}
TcMtextEntryWidgetE8 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE8)
tkgrid(TcMlabel1E8, TcMtextEntryWidgetE8)

TcMlabelTextE9 <- tclVar("Column of 2nd sequence position in 2nd file: ")
TcMlabel1E9 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE9))
if (length(con[which(con[,1]=="column of second position in second file"),2]) > 0){
	TcMtextEntryE9 <- tclVar(con[which(con[,1]=="column of second position in second file"),2])
}else{
	TcMtextEntryE9 <- tclVar("")
}
TcMtextEntryWidgetE9 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE9)
tkgrid(TcMlabel1E9, TcMtextEntryWidgetE9)

TcMlabelTextE0 <- tclVar("Column of p-values in 2nd file: ")
TcMlabel1E0 <- tklabel(graphicals2,text=tclvalue(TcMlabelTextE0))
if (length(con[which(con[,1]=="column of values in second file"),2]) > 0){
	TcMtextEntryE0 <- tclVar(con[which(con[,1]=="column of values in second file"),2])
}else{
	TcMtextEntryE0 <- tclVar("")
}
TcMtextEntryWidgetE0 <- tkentry(graphicals2,width=20,textvariable=TcMtextEntryE0)
tkgrid(TcMlabel1E0, TcMtextEntryWidgetE0)

TOK.but22 <- tkbutton(graphicals2,text="Start",command=function()calculate_co_mut_graphics_both(tbn5, txt, TcMtextEntryE, TcMtextEntryE1, TcMtextEntryE2, TcMtextEntryE3, TcMtextEntryE4, TcMtextEntryE5, TcMtextEntryE6, TcMtextEntryE7, TcMtextEntryE8, TcMtextEntryE9, TcMtextEntryE0))
tkgrid(TOK.but22)

#--------------------------------------------------------------------------------------------------------------------------------
tkgrid.rowconfigure(otbn3,0,weight=1, minsize=100)
tkgrid.columnconfigure(otbn3,0,weight=10, minsize=100)
tkgrid.columnconfigure(otbn3,1,weight=3, minsize=100)

tbn3 <- ttkframe(otbn3)
helptbn3 <- ttkframe(otbn3)
tkgrid(tbn3, helptbn3, sticky="snew")
tkgrid.columnconfigure(tbn3,0,weight=10)
tkgrid.columnconfigure(helptbn3,0,weight=1)

scr5 <- tkscrollbar(tbn3, command=function(...)tkyview(txt5,...),bg='white')
txt5 <- tktext(tbn3,height=1,width=5,bg='grey') 
xscr5 <- tkscrollbar(tbn3, command=function(...)tkxview(txt5,...),orient='horiz')
tkconfigure(txt5, yscrollcommand=function(...)tkset(scr5,...), xscrollcommand=function(...)tkset(xscr5,...), state="disabled")
tkpack(scr5,side='right',fill='y')
tkpack(txt5,expand='yes',fill='both')
tkpack(xscr5, side="bottom", fill="x")
sf3 <- tkcanvas(txt5, background="red")
tkwindow.create(txt5,'end', window=sf3)

qv <- tkframe(sf3,relief="groove",borderwidth=2)
tkgrid(qv, sticky="snew")

heading4 <- tklabel(qv,text="Calculate q-values from p-values", background="#FFCC99")
tkgrid(heading4, sticky="snew", columnspan=3)

button.widget_q <- tkbutton(qv,text="Select csv file",command=function()getcsvfile(qVlabelText1, 3, 1))

qVlabelText <- tclVar("Selected csv file: ")
qVlabel1 <- tklabel(qv,text=tclvalue(qVlabelText))
tkconfigure(qVlabel1,textvariable=qVlabelText)

qVlabelText1 <- tclVar("")
qVlabel2 <- tklabel(qv,text=tclvalue(qVlabelText1))
tkconfigure(qVlabel2,textvariable=qVlabelText1)
tkgrid(qVlabel1, qVlabel2, button.widget_q)

OK.but3 <- tkbutton(qv,text="Start",command=function()calculate_q_value(tbn3, txt))
tkgrid(OK.but3)
		#---------------------------------------------------------
epa <- tkframe(sf3,relief="groove",borderwidth=2)
tkgrid(epa, sticky="snew")

heading <- tklabel(epa,text="Possible features add on", background="#FFCC99")
tkgrid(heading, sticky="snew", columnspan=3)

button.widget_epa1 <- tkbutton(epa,text="Select FASTA file",command=function()getFASTAfile(pEFlabelText1, 4, 1))
pEFlabelText <- tclVar("Selected sequences: ")
pEFlabel1 <- tklabel(epa,text=tclvalue(pEFlabelText))
tkconfigure(pEFlabel1,textvariable=pEFlabelText)

pEFlabelText1 <- tclVar("")
pEFlabel2 <- tklabel(epa,text=tclvalue(pEFlabelText1))
tkconfigure(pEFlabel2,textvariable=pEFlabelText1)
tkgrid(pEFlabel1, pEFlabel2, button.widget_epa1)

button.widget_epa2 <- tkbutton(epa,text="Select csv file",command=function()getcsvfile(pEFlabelText1A, 4, 2))
pEFlabelTextA <- tclVar("Selected result from point mutations: ")
pEFlabel1A <- tklabel(epa,text=tclvalue(pEFlabelTextA))
tkconfigure(pEFlabel1A,textvariable=pEFlabelTextA)

pEFlabelText1A <- tclVar("")
pEFlabel2A <- tklabel(epa,text=tclvalue(pEFlabelText1A))
tkconfigure(pEFlabel2A,textvariable=pEFlabelText1A)
tkgrid(pEFlabel1A, pEFlabel2A, button.widget_epa2)

button.widget_epa3 <- tkbutton(epa,text="Select csv file",command=function()getcsvfile(pEFlabelText1B, 4, 3))
pEFlabelTextB <- tclVar("Selected results from tuple mutations: ")
pEFlabel1B <- tklabel(epa,text=tclvalue(pEFlabelTextB))
tkconfigure(pEFlabel1B,textvariable=pEFlabelTextB)

pEFlabelText1B <- tclVar("")
pEFlabel2B <- tklabel(epa,text=tclvalue(pEFlabelText1B))
tkconfigure(pEFlabel2B,textvariable=pEFlabelText1B)
tkgrid(pEFlabel1B, pEFlabel2B, button.widget_epa3)

pElabelText2 <- tclVar("p-value: ")
pElabel2G <- tklabel(epa,text=tclvalue(pElabelText2))
if (length(con[which(con[,1]=="p-value addon"),2]) > 0){
	textEntry_p.value <- tclVar(con[which(con[,1]=="p-value addon"),2])
}else{
	textEntry_p.value <- tclVar("")
}
textEntry_Wp.value <- tkentry(epa,width=20,textvariable=textEntry_p.value)
tkgrid(pElabel2G, textEntry_Wp.value)

OK.bute2 <- tkbutton(epa,text="Start",command=function()calculate_pos_epi_further(tbn3, textEntry_p.value, txt))
tkgrid(OK.bute2)
		#-----------------------------------------------------
#Founder elem
fe <- tkframe(sf3,relief="groove",borderwidth=2)
tkgrid(fe, sticky="snew")

heading3 <- tklabel(fe,text="Founder effect eliminator", background="#FFCC99")
tkgrid(heading3, sticky="snew", columnspan=3)

button.widget_fe1 <- tkbutton(fe,text="Select FASTA file",command=function()getFASTAfile(felabelText1, 6, 1))
felabelText <- tclVar("Selected sequences: ")
felabel1 <- tklabel(fe,text=tclvalue(felabelText))
tkconfigure(felabel1,textvariable=felabelText)
felabelText1 <- tclVar("")
felabel2 <- tklabel(fe,text=tclvalue(felabelText1))
tkconfigure(felabel2,textvariable=felabelText1)
tkgrid(felabel1, felabel2, button.widget_fe1)

button.widget_fe2 <- tkbutton(fe,text="Select csv file",command=function()getcsvfile(felabelText12, 6, 2))
felabelText2 <- tclVar("Selected results from tuple mutations: ")
felabel12 <- tklabel(fe,text=tclvalue(felabelText2))
tkconfigure(felabel12,textvariable=felabelText2)
felabelText12 <- tclVar("")
felabel22 <- tklabel(fe,text=tclvalue(felabelText12))
tkconfigure(felabel22,textvariable=felabelText12)
tkgrid(felabel12, felabel22, button.widget_fe2)

button.widget_fe3 <- tkbutton(fe,text="Select nexus file",command=function()getnexusfile(felabelText13, 6, 3))
felabelText3 <- tclVar("Selected tree of loaded sequences: ")
felabel13 <- tklabel(fe,text=tclvalue(felabelText3))
tkconfigure(felabel13,textvariable=felabelText3)
felabelText13 <- tclVar("")
felabel23 <- tklabel(fe,text=tclvalue(felabelText13))
tkconfigure(felabel23,textvariable=felabelText13)
tkgrid(felabel13, felabel23, button.widget_fe3)

felabelTextE <- tclVar("Ratio in subtree: ")
felabel1E <- tklabel(fe,text=tclvalue(felabelTextE))
if (length(con[which(con[,1]=="Ratio in subtree"),2]) > 0){
	fetextEntryE <- tclVar(con[which(con[,1]=="Ratio in subtree"),2])
}else{
	fetextEntryE <- tclVar("")
}
fetextEntryWidgetE <- tkentry(fe,width=20,textvariable=fetextEntryE)
tkgrid(felabel1E, fetextEntryWidgetE)

OK.but62 <- tkbutton(fe,text="Start",command=function()calculate_founder(tbn3, fetextEntryE, txt))
tkgrid(OK.but62)

tkgrid.configure(button.widget_co1, column=7)
tkgrid.configure(button.widget_co2, column=7)
tkgrid.configure(cMlabel2, columnspan=5)
tkgrid.configure(cMlabel2A, columnspan=5)
tkgrid.configure(cMtextEntryWidgetC, columnspan=5)
tkgrid.configure(cMtextEntryWidgetD, columnspan=5)
tkgrid.configure(cbm, columnspan=5)
#### Rewrite shared mutations----------------------------
rsm <- tkframe(sf3,relief="groove",borderwidth=2)
tkgrid(rsm, sticky="snew")

heading3 <- tklabel(rsm,text="Rewrite tuple analysis result for tartan", background="#FFCC99")
tkgrid(heading3, sticky="snew", columnspan=7)

button.widget_rsm2 <- tkbutton(rsm,text="Select csv file",command=function()getcsvfile(rsmlabelText12, 10, 1))
rsmlabelText2 <- tclVar("Selected results from tuple analysis: ")
rsmlabel12 <- tklabel(rsm,text=tclvalue(rsmlabelText2))
tkconfigure(rsmlabel12,textvariable=rsmlabelText2)
rsmlabelText12 <- tclVar("")
rsmlabel22 <- tklabel(rsm,text=tclvalue(rsmlabelText12))
tkconfigure(rsmlabel22,textvariable=rsmlabelText12)
tkgrid(rsmlabel12, rsmlabel22, button.widget_rsm2)

rsmlabelTextE5 <- tclVar("Column of 1st sequence position: ")
rsmlabel1E5 <- tklabel(rsm,text=tclvalue(rsmlabelTextE5))
if (length(con[which(con[,1]=="rw column of first position"),2]) > 0){
	rsmtextEntryE5 <- tclVar(con[which(con[,1]=="rw column of first position"),2])
}else{
	rsmtextEntryE5 <- tclVar("")
}
rsmtextEntryWidgetE5 <- tkentry(rsm,width=20,textvariable=rsmtextEntryE5)
tkgrid(rsmlabel1E5, rsmtextEntryWidgetE5)

rsmlabelTextE6 <- tclVar("Column of 2nd sequence position: ")
rsmlabel1E6 <- tklabel(rsm,text=tclvalue(rsmlabelTextE6))
if (length(con[which(con[,1]=="rw column of second position"),2]) > 0){
	rsmtextEntryE6 <- tclVar(con[which(con[,1]=="rw column of second position"),2])
}else{
	rsmtextEntryE6 <- tclVar("")
}
rsmtextEntryWidgetE6 <- tkentry(rsm,width=20,textvariable=rsmtextEntryE6)
tkgrid(rsmlabel1E6, rsmtextEntryWidgetE6)

rsmlabelTextE7 <- tclVar("Column of p-values: ")
rsmlabel1E7 <- tklabel(rsm,text=tclvalue(rsmlabelTextE7))
if (length(con[which(con[,1]=="rw column of values"),2]) > 0){
	rsmtextEntryE7 <- tclVar(con[which(con[,1]=="rw column of values"),2])
}else{
	rsmtextEntryE7 <- tclVar("")
}
rsmtextEntryWidgetE7 <- tkentry(rsm,width=20,textvariable=rsmtextEntryE7)
tkgrid(rsmlabel1E7, rsmtextEntryWidgetE7)

rsmlabelTextA <- tclVar("Csv seperator: ")
rsmlabel1A <- tklabel(rsm,text=tclvalue(rsmlabelTextA))
if (length(con[which(con[,1]=="rw csv seperator"),2]) > 0){
	rsmtextEntryA <- tclVar(con[which(con[,1]=="rw csv seperator"),2])
}else{
	rsmtextEntryA <- tclVar("")
}
rsmtextEntryWidgetA <- tkentry(rsm,width=5,textvariable=rsmtextEntryA)
tkgrid(rsmlabel1A, rsmtextEntryWidgetA)

rsmlabelTextB <- tclVar("Threshold: ")
rsmlabel1B <- tklabel(rsm,text=tclvalue(rsmlabelTextB))
if (length(con[which(con[,1]=="rw threshold"),2]) > 0){
	rsmtextEntryB <- tclVar(con[which(con[,1]=="rw threshold"),2])
}else{
	rsmtextEntryB <- tclVar("")
}
rsmtextEntryWidgetB <- tkentry(rsm,width=20,textvariable=rsmtextEntryB)
tkgrid(rsmlabel1B, rsmtextEntryWidgetB)

OK.but72 <- tkbutton(rsm,text="Start",command=function()rewrite_shared_mutations(tbn3, rsmtextEntryE5, rsmtextEntryE6, rsmtextEntryE7, rsmtextEntryA, rsmtextEntryB, txt))
tkgrid(OK.but72)

#--------------------------------------------------------HELP-TEXT-----------------------
txthelptbn <- tktext(helptbn,bg="white", font="courier", width=20, wrap="word")
tkgrid(txthelptbn)

tkgrid.rowconfigure(txthelptbn,0,weight=1)
tkgrid.columnconfigure(txthelptbn,0,weight=1)
tkgrid.configure(txthelptbn,sticky='nswe')

helptbntext <- paste("1. Make sure your data is in FASTA format and your feature is properly formated:\n",
"1.1. If you have HLA types as features (HLA-A and HLA-B) all types have to be listed in the comment line of each sequence on the same position.\n",
"1.2. If you have just one feature: Please separate this feature from the rest of the comment line with an \";\"\n",
"2. It is important for the orPlot function, that you enter the correct columns of everything. For this you have to count them inside your csv result file from the first block."
)

tryCatch({tkinsert(txthelptbn,"end",helptbntext)}, error = function(err) {})
tkconfigure(txthelptbn, state="disabled")
#-----
txthelptbn2 <- tktext(helptbn2,bg="white", font="courier", width=20, wrap="word")
tkgrid(txthelptbn2)

tkgrid.rowconfigure(txthelptbn2,0,weight=1)
tkgrid.columnconfigure(txthelptbn2,0,weight=1)
tkgrid.configure(txthelptbn2,sticky='nswe')

helptbn2text <- paste("1.1 Find those position pairs above the chosen threshold which have an p-value above the chosen threshold.\n",
"You can choose whether to use feature specific analysis or not\n",
"1.2 Find those position tuple(more than two members possible) which have an p-value above the chosen threshold from the analysis of point mutations and features."

)

tryCatch({tkinsert(txthelptbn2,"end",helptbn2text)}, error = function(err) {})
tkconfigure(txthelptbn2, state="disabled")
#-----
txthelptbn5 <- tktext(helptbn5,bg="white", font="courier", width=20, wrap="word")
tkgrid(txthelptbn5)

tkgrid.rowconfigure(txthelptbn5,0,weight=1)
tkgrid.columnconfigure(txthelptbn5,0,weight=1)
tkgrid.configure(txthelptbn5,sticky='nswe')

helptbn5text <- paste("1. For the tartan plot you can compare the results of two different pair mutation files. Be aware that for the additional ticks and notes in the plot you have to enter values which are inside the range of the sequence length which created this data.\n"
)

tryCatch({tkinsert(txthelptbn5,"end",helptbn5text)}, error = function(err) {})
tkconfigure(txthelptbn5, state="disabled")
#-----
txthelptbn3 <- tktext(helptbn3,bg="white", font="courier", width=20, wrap="word")
tkgrid(txthelptbn3)

tkgrid.rowconfigure(txthelptbn3,0,weight=1)
tkgrid.columnconfigure(txthelptbn3,0,weight=1)
tkgrid.configure(txthelptbn3,sticky='nswe')

helptbn3text <- paste("1. The file for added q-value column has to have a column with \"p_value\" (not p-value!!) as heading.\n",
"2. For analysis of a possible founders effect DON'T use an pair mutation result which differentiates between identifiers.\n" 

)

tryCatch({tkinsert(txthelptbn3,"end",helptbn3text)}, error = function(err) {})
tkconfigure(txthelptbn3, state="disabled")
#--------------------------------------------------------TOOL-TIPS-----------------------
tk2tip(cb01, "If checked FASTA file has to have nucleotides, if not amino acids.")
tk2tip(icb01, "If checked data has only one feature (like tropism), if not two (like HLA types).")

tk2tip(button.widget_epi1, "Loads FASTA file")
tk2tip(button.widget_epi2, "Loads file with known epitopes. See example files.")
tk2tip(button.widget_epi3, "Loads file with known binding motifs. See example files.")
tk2tip(button.widget_epi4, "Load reference sequence for own position corrections, if you delete position columns due to lots of gaps.")
tk2tip(textEntryWidgetC, "The minimum number of members which have the same feature to be included in the analysis.")
tk2tip(textEntryWidgetD, "The optical level height of \"significant\" results in the graphical output.")
tk2tip(textEntryWidgetE, "The optical level height of \"star level\" results in the graphical output. Here he corrected column is used.")
tk2tip(textEntryWidgetF, "p-value below which a position is tested of binding motifs (yellow line).")
tk2tip(textEntryWidgetG, "p-value below which a position is counted for sliding window (black line).")
tk2tip(OK.bute1, "Starts the calculation.")

tk2tip(textEntryWidgetH, "Position of feature-A_first start in FASTA comment line.")
tk2tip(textEntryWidgetH1, "Position of feature-A_first end in FASTA comment line.")
tk2tip(textEntryWidgetH2, "Position of feature-A_second start in FASTA comment line.")
tk2tip(textEntryWidgetH12,"Position of feature-A_second end in FASTA comment line.")
tk2tip(textEntryWidgetHB, "Position of feature-B_first start in FASTA comment line.")
tk2tip(textEntryWidgetH1B, "Position of feature-B_first end in FASTA comment line.")
tk2tip(textEntryWidgetH2B, "Position of feature-B_second start in FASTA comment line.")
tk2tip(textEntryWidgetH12B, "Position of feature-B_second end in FASTA comment line.")
tk2tip(cb0e, "Should a phylogenetic comparison be made?")
tk2tip(textEntryWidgetKIB, "Which matrix should be applied for phylogenetic comparison?")
tk2tip(ddb, "Input can be: \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\", \"none\".")
tk2tip(textEntryWidgetH42B, "The size of the sliding window for graphical output. Must be greater than 1.")

tk2tip(Bbutton.widget_epi, "Loads FASTA file.")
tk2tip(BtextEntryWidgetA, "The separator used in the csv file.")
tk2tip(BtextEntryWidgetB, "The number of cases aka different identifiers for your feature(s) in the csv file.")
tk2tip(BtextEntryWidgetC, "The first column with an odds.ratio.")
tk2tip(BtextEntryWidgetD, "The first column with p-values.")
tk2tip(BtextEntryWidgetE, "The column with the name for the case in the first row.")
tk2tip(BtextEntryWidgetF, "The frequency of columns aka the number of columns for one identifier.")
tk2tip(BtextEntryWidgetG, "The first column of amino acids (if existed).")
tk2tip(BtextEntryWidgetH, "The maximum value on the y-axis (on both sides).")
tk2tip(BtextEntryWidgetI, "The main interval on the y-axis.")
tk2tip(BtextEntryWidgetJ, "If a color should be added, enter column with values between 0 and 1.")
tk2tip(BOK.bute2, "Starts the plot making.")

tk2tip(button.widget_epa1, "Loads FASTA file.")
tk2tip(button.widget_epa2, "Loads file with point mutation results.")
tk2tip(button.widget_epa3, "Loads file with co-mutation results.")
tk2tip(textEntry_Wp.value, "p-value for selecting results from co-mutation.")
tk2tip(OK.bute2, "Starts the calculation.")
#----
tk2tip(cb1, "Should analysis include different allels or compare to consensus?")
#----
tk2tip(cb, "Should analysis include different allels or compare to consensus?")
tk2tip(button.widget_co1, "Loads FASTA file.")
tk2tip(button.widget_co2, "Loads file with Consensus of sequences in FASTA file.")
tk2tip(cMtextEntryWidgetC, "The minimum number of patients which have the same feature type to be included in the analysis.")
tk2tip(cMtextEntryWidgetD, "p-value below possible co-mutated positions are included in the analysis.")
tk2tip(ddbc, "Input can be: \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\", \"none\".")
tk2tip(cbm, "Should calculation be made with more than one core?")
tk2tip(OK.but21, "Starts the calculation.")

tk2tip(textEntryWidgetHC, "Position of feature-A_first start in sequence description.")
tk2tip(textEntryWidgetH1C, "Position of feature-A_first end in sequence description.")
tk2tip(textEntryWidgetH2C, "Position of feature-A_second start in sequence description.")
tk2tip(textEntryWidgetH12C,"Position of feature-A_second end in sequence description.")
tk2tip(textEntryWidgetHBC, "Position of feature-B_first start in sequence description.")
tk2tip(textEntryWidgetH1BC, "Position of feature-B_first end in sequence description.")
tk2tip(textEntryWidgetH2BC, "Position of feature-B_second start in sequence description.")
tk2tip(textEntryWidgetH12BC, "Position of feature-B_second end in sequence description.")

tk2tip(button.widget_com, "Loads co-mutation results. Be careful with checking \"allels\" or \"no allels\".")
tk2tip(cMtextEntryWidgetE, "p-value for optical differentiation.")
tk2tip(OK.but22, "Starts the calculation.")
#----
tk2tip(Tbutton.widget_com, "The first file with two position columns and p-values for those.")
tk2tip(Tbutton.widget_com2, "The second file with two position columns and p-values for those.")
tk2tip(Tbutton.widget_com3, "A file with distance matrix values.")
tk2tip(TcMtextEntryWidgetE, "The optical space in px between the blocks.")
tk2tip(TcMtextEntryWidgetE1, "The colors of the plot, ranging from lowest to highest. Given as list: \"color1, color2,")
tk2tip(TcMtextEntryWidgetE2, "The positions for names for interesting positions. Given as list: \"pos1, pos2,..\"")
tk2tip(TcMtextEntryWidgetE3, "The names for interesting positions. Given as list: \"name1, name2,..\"")
tk2tip(TcMtextEntryWidgetE4, "The ticks on x and y axis.")
tk2tip(TcMtextEntryWidgetE5, "The column of the first position in the first file.")
tk2tip(TcMtextEntryWidgetE6, "The column of the second position in the first file.")
tk2tip(TcMtextEntryWidgetE7, "The column of the value in the second file.")
tk2tip(TcMtextEntryWidgetE8, "The column of the first position in the second file.")
tk2tip(TcMtextEntryWidgetE9, "The column of the second position in the second file.")
tk2tip(TcMtextEntryWidgetE0, "The column of the value in the first file.")
tk2tip(TOK.but22, "Starts the plotmaking")
#----
tk2tip(button.widget_sm1, "The fasta file with sequences")
tk2tip(button.widget_sm2, "The result file from point mutation analysis or other file with positions and (p-)values")
tk2tip(smtextEntryWidgetA, "The threshold below which a position should be taken for calculation.")
tk2tip(smtextEntryWidgetB, "The minimal number of tuple for position tuple.")
tk2tip(smtextEntryWidgetC, "The maximal number of tuple for position tuple.")
tk2tip(smtextEntryWidgetD, "The column in which the identifier is column name.")
tk2tip(smtextEntryWidgetE, "The column of p-values for this identifier.")
tk2tip(smtextEntryWidgetF, "The column of amino acids/nucleotides for this identifier.")
tk2tip(smtextEntryWidgetH, "Position of feature-A_first start in sequence description.")
tk2tip(smtextEntryWidgetH1, "Position of feature-A_first end in sequence description.")
tk2tip(smtextEntryWidgetH2, "Position of feature-A_second start in sequence description.")
tk2tip(smtextEntryWidgetH12,"Position of feature-A_second end in sequence description.")
tk2tip(smtextEntryWidgetHB, "Position of feature-B_first start in sequence description.")
tk2tip(smtextEntryWidgetH1B, "Position of feature-B_first end in sequence description.")
tk2tip(smtextEntryWidgetH2B, "Position of feature-B_second start in sequence description.")
tk2tip(smtextEntryWidgetH12B, "Position of feature-B_second end in sequence description.")
#tk2tip(as.character(tclvalue(icbValue01)), "")
tk2tip(smtextEntryWidgetX , "The identifier which should be analyzed (if there is a single one).")
tk2tip(OK.but4, "Starts the calculation.")
#----
tk2tip(button.widget_q, "Load file with \"p_value\" column.")
tk2tip(OK.but3, "Starts the calculation.")

tk2tip(button.widget_fe1, "Loads FASTA file")
tk2tip(button.widget_fe2, "Loads co-mutation results. Not possible in combination with \"allels\".")
tk2tip(button.widget_fe3, "Loads tree file in nexus format.")
tk2tip(fetextEntryWidgetE, "Rate of co-mutation is in subtree vs is not in subtree.")
tk2tip(OK.but62, "Starts the calculation")

#---------------------------------------------------------
tkbind(tt, "<Destroy>", function() tkdestroy(tn))
tkselect(tn, 0)
#-------------FOR CONFIG--------------------------------------------
z <- list()
z$cbValue01 <- list(cbValue01, "DNA")
z$cbValue_0e <- list(cbValue_0e, "Phylogenetic comparison")
z$textEntryC <- list(textEntryC, "Number of patients threshold")
z$textEntryD <- list(textEntryD, "Height of horizontal bar")
z$textEntryE <- list(textEntryE, "Height of star level")
z$textEntryF <- list(textEntryF, "Level for possible epitope")
z$textEntryG <- list(textEntryG, "Window Threshold")
z$textEntryH <- list(textEntryH, "HLA-A1")
z$textEntryH1 <- list(textEntryH1, "HLA-A2")
z$textEntryH2 <- list(textEntryH2, "HLA-A3")
z$textEntryH12 <- list(textEntryH12, "HLA-A4")
z$textEntryHB <- list(textEntryHB, "HLA-B1")
z$textEntryH1B <- list(textEntryH1B, "HLA-B2")
z$textEntryH2B <- list(textEntryH2B, "HLA-B3")
z$textEntryH12B <- list(textEntryH12B, "HLA-B4")
z$textEntryIB <- list(textEntryIB, "p-value correction")
z$textEntryKIB <- list(textEntryKIB, "Matrix for phylo")
z$icbValue01 <- list(icbValue01, "One identifier")
z$textEntryH42B <- list(textEntryH42B, "window_size")
z$textEntry_p.value <- list(textEntry_p.value, "p-value addon")
z$BtextEntryA <- list(BtextEntryA, "csv seperator")
z$BtextEntryB <- list(BtextEntryB, "number of cases")
z$BtextEntryC <- list(BtextEntryC, "position of odds.ratio")
z$BtextEntryD <- list(BtextEntryD, "position of p.values")
z$BtextEntryE <- list(BtextEntryE, "name for the y-axis label")
z$BtextEntryF <- list(BtextEntryF, "frequency")
z$BtextEntryG <- list(BtextEntryG, "position of amino acid")
z$BtextEntryH <- list(BtextEntryH, "maximum value of y-axis")
z$BtextEntryI <- list(BtextEntryI, "intervall of y-axis")
z$BtextEntryJ <- list(BtextEntryJ, "add color information")
z$cMtextEntryC <- list(cMtextEntryC, "Number of patients threshold_cm")
z$cMtextEntryD <- list(cMtextEntryD, "level for significance")
z$cmtextEntryIB <- list(cmtextEntryIB, "cmp-value correction")
z$cbValuem <- list(cbValuem, "More than one core")
z$textEntryHC <- list(textEntryHC, "cmHLA-A1")
z$textEntryH1C <- list(textEntryH1C, "cmHLA-A2")
z$textEntryH2C <- list(textEntryH2C, "cmHLA-A3")
z$textEntryH12C <- list(textEntryH12C,"cmHLA-A4")
z$textEntryHBC <- list(textEntryHBC, "cmHLA-B1")
z$textEntryH1BC <- list(textEntryH1BC, "cmHLA-B2")
z$textEntryH2BC <- list(textEntryH2BC, "cmHLA-B3")
z$textEntryH12BC <- list(textEntryH12BC, "cmHLA-B4")
z$fetextEntryE <- list(fetextEntryE, "Ratio in subtree")
z$cMtextEntryE <- list(cMtextEntryE, "cMlevel for significance")
z$TcMtextEntryE <- list(TcMtextEntryE, "space between blocks")
z$TcMtextEntryE1 <- list(TcMtextEntryE1, "colors of the plot")
z$TcMtextEntryE2 <- list(TcMtextEntryE2, "name positions")
z$TcMtextEntryE3 <- list(TcMtextEntryE3, "names")
z$TcMtextEntryE4 <- list(TcMtextEntryE4, "ticks")
z$TcMtextEntryE5 <- list(TcMtextEntryE5, "column of first position in first file")
z$TcMtextEntryE6 <- list(TcMtextEntryE6, "column of second position in first file")
z$TcMtextEntryE7 <- list(TcMtextEntryE7, "column of values in first file")
z$TcMtextEntryE8 <- list(TcMtextEntryE8, "column of first position in second file")
z$TcMtextEntryE9 <- list(TcMtextEntryE9, "column of second position in second file")
z$TcMtextEntryE0 <- list(TcMtextEntryE0, "column of values in second file")
z$smtextEntryA <- list(smtextEntryA, "threshold")
z$smtextEntryB <- list(smtextEntryB, "min_number_of_ele_in_tupel")
z$smtextEntryC <- list(smtextEntryC, "max_number_of_ele_in_tupel")
z$smtextEntryD <- list(smtextEntryD, "column")
z$smtextEntryE <- list(smtextEntryE, "column of values")
z$smtextEntryF <- list(smtextEntryF, "column of aas")
z$smtextEntryH <- list(smtextEntryH, "smHLA-A1")
z$smtextEntryH1 <- list(smtextEntryH1, "smHLA-A2")
z$smtextEntryH2 <- list(smtextEntryH2, "smHLA-A3")
z$smtextEntryH12 <- list(smtextEntryH12,"smHLA-A4")
z$smtextEntryHB <- list(smtextEntryHB, "smHLA-B1")
z$smtextEntryH1B <- list(smtextEntryH1B, "smHLA-B2")
z$smtextEntryH2B <- list(smtextEntryH2B, "smHLA-B3")
z$smtextEntryH12B <- list(smtextEntryH12B, "smHLA-B4")
z$smtextEntryidet <- list(smtextEntryidet, "Identifier")
z$rsmtextEntryE5 <- list(rsmtextEntryE5, "rw column of first position")
z$rsmtextEntryE6 <- list(rsmtextEntryE6, "rw column of second position")
z$rsmtextEntryE7 <- list(rsmtextEntryE7, "rw column of values")
z$rsmtextEntryA <- list(rsmtextEntryA, "rw csv seperator")
z$rsmtextEntryB <- list(rsmtextEntryB, "rw threshold")
.GlobalEnv[["z"]] <- z
#---------------------------------------------------------
},ex=function(){
	SeqFeatR_GUI()
})

#SeqFeatR_GUI()
