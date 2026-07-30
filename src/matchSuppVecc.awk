#! /bin/awk -f

#expects names, colindex

BEGIN{
    chr=1;
    pos=2;
    tempcol=3;
    svid=4;
    ref=5;
    alt=6;
    svlen=7;
    suppvecc=colindex;
    svtype=9;
    OFS="\t"
    
    split(names, nameArray, ",");
}
{
   # Get the binary string from column 8
    binary_string = $suppvecc;
    
    # Check if the binary string is all zeros
    if (binary_string == "" || binary_string == "." || binary_string ~ /^\.?0+$/) {
        names_string = ".";
    } else {
        # Initialize the names string
        names_string = "";
        
        # Loop over each character in the binary string
        for (i = 1; i <= length(binary_string); i++) {
            if (substr(binary_string, i, 1) == "1") {
                if (names_string == "") {
                    names_string = nameArray[i];
                } else {
                    names_string = names_string "_" nameArray[i];
                }
            }
        }
    }
    
    # Replace column 8 with the names string
    $suppvecc = names_string;
    
    # Print the modified row
    print $0;
}
