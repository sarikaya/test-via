requirements = """
For the library column, please use the entered names while creating a collection.
"""

def quo2Lst (str){
    list = str.replaceAll("'","").split(",")
    return list
}
sampleIndex = -1
groupIndex = -1
groupSep = ","
errorText = ""
fileNamesInPairFile = []
fileNamesInGroupFile = []
if (params.metadata && params.metadata != "NA"){
    initialrunFileNames = new File("./../initialrun/.file_name")
    fileNames =[]
    if (initialrunFileNames.exists()){
        file_name = initialrunFileNames.text
        fileNames = quo2Lst(file_name)
    }
    
    
    println(params.metadata)
    new FileReader(params.metadata).eachLine { line, number ->
        

        if (number == 1 && !line.contains("library")){
            errorText += "$line \n"
            errorText += "ERROR: library column is missing in metadata file\n"
        }
        if (number == 1 && !line.contains("group")){
            errorText += "$line \n"
            errorText += "ERROR: library column is missing in the metadata file\n"
        }
        if (number == 1 && line.contains(" ")){
            errorText += "$line \n"
            errorText += "ERROR: Space is not allowed for the column names in the metadata file \n"
        }
        if (number == 1 && line.contains("-")){
            errorText += "$line \n"
            errorText += "ERROR: Dash/Hyphen is not allowed for the column names in the metadata file \n"
        }
        
        groupVals = line.split(groupSep)
        if (number == 1) sampleIndex = groupVals.findIndexOf { it.endsWith("library")}
        if (sampleIndex >-1 && number != 1 && sampleIndex < groupVals.size()) fileNamesInGroupFile.push(groupVals[sampleIndex])
    }
if (params.metadata_paired && params.metadata_paired != "NA"){
    
    new FileReader(params.metadata_paired).eachLine { line, number ->
        
        if (line.contains("\t")){
          groupSep = "\t" 
        }
        if (number == 1 && !line.contains("p3_library")){
            errorText += "$line \n"
            errorText += "ERROR: p3_library column is missing in metadata_paired file\n"
        }
        if (number == 1 && !line.contains("p5_library")){
            errorText += "$line \n"
            errorText += "ERROR: p5_library column is missing in the metadata_paired file\n"
        }
        
        
        groupVals = line.split(groupSep)
        if (number == 1) sampleIndex3 = groupVals.findIndexOf { it ==~ /p3_library/ }
        if (number == 1) sampleIndex5 = groupVals.findIndexOf { it ==~ /p5_library/ }
        if (sampleIndex3 >-1 && number != 1 && sampleIndex3 < groupVals.size()) fileNamesInPairFile.push(groupVals[sampleIndex3])
        if (sampleIndex5 >-1 && number != 1 && sampleIndex5 < groupVals.size()) fileNamesInPairFile.push(groupVals[sampleIndex5])

    }
    }
    
    difference = []
    if (fileNames.size() > 0){
        fileNamesInGroupFile.eachWithIndex{ it, i -> {
            it2=it
            ind = fileNames.findIndexOf { it == it2 }
            if (ind < 0 && it2 != "") {
                errorText += "ERROR: library \"$it2\" is not found in the collection.\n"
                difference.push(it2)
            }
        }
        }
         fileNamesInPairFile.eachWithIndex{ it, i -> {
            it2=it
            ind = fileNames.findIndexOf { it == it2 }
            if (ind < 0 && it2 != "") {
                errorText += "ERROR: paired library \"$it2\" is not found in the collection.\n"
                difference.push(it2)
            }
        } 
            
    }
    }

    if (difference.size() > 0){
        errorText += "Collection Sample Names: " +fileNames + "\n"
        errorText += "Metadata File Sample Names: " +fileNamesInGroupFile + "\n"
        errorText += "Paired Metadata File Sample Names: " +fileNamesInPairFile + "\n"
    }
    
    
    
    if (errorText != ""){
        println(errorText)
        System.exit(1)
    }
    
}

chromosome_list = params.Check_and_Build_Module_Add_chromosomes_to_genome_gtf.chromosome_list
initialList = chromosome_list.split("\\s+")
warnlist = []
for (int i = 0; i < initialList.size(); i++) {
    if (initialList[i] ==~ /SQ-\d{6}$/){
    } else {
        warnlist.push(initialList[i])
    }
}

if (chromosome_list?.trim() && !warnlist.isEmpty()){
    println ("ERROR: chromosome_list should be space-separated SQ IDs which must be in the format: SQ-000001 SQ-000002. Following item(s) are not matching the expected format")
    println (warnlist)
    System.exit(1)
}

