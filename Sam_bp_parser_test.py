str1 = "HISEQ05:409:C5VPLACXX:5:1210:12718:20738        161     JMPD004_index8_Contig1  1       70      74M     =       310     402     ATTTTCAGTCCCAGATTATGCCATGATCTAGTATCTGCCAATTTGGACAGGGGCATTAGTTCACAATTTCCCTT      =;?DDDDBAA>DD?C<?B:<+<<3<<33C?E994CFIAEI*:?D<B3B3?=BDDAD@BD==4BC@88@=777A@      MC:Z:93M        MD:Z:74 NM:i:0  MQ:i:70 UQ:i:1  AS:i:1  RG:Z:JMPD004_index8     PG:Z:novoalign"
str2 = "HISEQ05:409:C5VPLACXX:5:1308:12865:94980        99      JMPD004_index8_Contig1  1       70      99M1H   =       173     272     ATTTTCAGTCCCAGATTATGCCATGATCTAGTATCTGCCAATTTGGACAGGGGCATTAGTTCACAATTTCCCTTGTTGTTGCCAACCTTTTAATCGGCA     CCCFFFFFHHHHHJJJJJJJJJJJJJJIJJJJJJJJJJJIGIJJJIJIIJJJJIJJIJIIIJIJJJJJJGGIIJJHHHFHFFFFFDDEEDDDDDD@BDD     MC:Z:100M       MD:Z:99 AM:i:70 NM:i:0  SM:i:70 MQ:i:70 PQ:i:29 UQ:i:0  AS:i:0  RG:Z:JMPD004_index8     PG:Z:novoalign"
str3 = "HISEQ05:409:C5VPLACXX:5:1311:14186:72134        163     JMPD004_index8_Contig1  1       70      2S169M1H   =       178     277     ATTTTCAGTCCCAGATTATGCCATGATCTAGTATCTGCCAATTTGGACAGGGGCATTAGTTCACAATTTCCCTTGTTGTTGCCAACCTTTTATTCGGCA     @CCFFFDFHHHGHDGIJEIIJJIIJJIIJIJHEGIIIJJJJIGIJJJJJIJIIJJJJIGHIIIJJJJJJJIJJJJHHFHHFFFFFDDEEDDDDDEDBDD     MC:Z:100M       MD:Z:92A6       AM:i:70 NM:i:1  SM:i:70 MQ:i:70 PQ:i:63 UQ:i:30 AS:i:30 RG:Z:JMPD004_index8     PG:Z:novoalign"
str4 = "HISEQ05:409:C5VPLACXX:5:1213:10173:28739        163     JMPD004_index8_Contig1  4       70      2S2M8I756M1H   =       199     295     TTCAGTCCCAGATTATGCCATGATCTAGTATCTGCCAATTTGGACAGGGGCATTAGTTCACAATTTCCCTTGTTGTTGCCAACCTTTTAATCGGCAACA     C@CFFFFFHHH;DGIHGIJJJJIIIJJJHIIGIJIJJJGHIHJDIIGIJJGIIIJJIJJJJJJJJIJHIIJJIIHHHHGFFFFDEEDEEDDDDDDDDDD     MC:Z:100M       MD:Z:99 AM:i:70 NM:i:0  SM:i:70 MQ:i:70 PQ:i:50 UQ:i:0  AS:i:0  RG:Z:JMPD004_index8     PG:Z:novoalign"
str5 = "HISEQ05:409:C5VPLACXX:5:1213:10173:28739        163     JMPD004_index8_Contig1  4       70      2M8I26M1H2M   =       199     295     TTCAGTCCCAGATTATGCCATGATCTAGTATCTGCCAATTTGGACAGGGGCATTAGTTCACAATTTCCCTTGTTGTTGCCAACCTTTTAATCGGCAACA     C@CFFFFFHHH;DGIHGIJJJJIIIJJJHIIGIJIJJJGHIHJDIIGIJJGIIIJJIJJJJJJJJIJHIIJJIIHHHHGFFFFDEEDEEDDDDDDDDDD     MC:Z:100M       MD:Z:99 AM:i:70 NM:i:0  SM:i:70 MQ:i:70 PQ:i:50 UQ:i:0  AS:i:0  RG:Z:JMPD004_index8     PG:Z:novoalign"
list1 = []

list1.append(str1)
list1.append(str2)
list1.append(str3)
list1.append(str4)
list1.append(str5)


def base_counter(x):
    i = int(0)
    line_bp = int()
    for char in x:
        if i <= int(1):
            if char == 'M':
                h = i - int(1)
                if base_str[h].isdigit():
                    bp = int(base_str[h])
                    #print bp
                    line_bp += bp
                    
        elif i <= int(2) and i >= int(1):
            if char == 'M':
                h = i - int(1)
                g = i - int(2)
                if base_str[h].isdigit() and base_str[g].isdigit():
                    bp = int( str(base_str[g]) + str(base_str[h]))
                    #print bp
                    line_bp += bp
                elif base_str[h].isdigit() and base_str[g].isdigit() == False:
                    bp = int(base_str[h])
                    #print bp
                    line_bp += bp
                    
        elif i > int(2):
            if char == 'M':
                h = i - int(1)
                g = i - int(2)
                f = i - int(3)
                if base_str[h].isdigit() and base_str[g].isdigit() and base_str[f].isdigit():
                    bp = int( str(base_str[f]) + str(base_str[g]) + str(base_str[h]))
                    #print bp
                    line_bp += bp
                elif base_str[h].isdigit() and base_str[g].isdigit() and base_str[f].isdigit() == False:
                    bp = int(str(base_str[g]) + str(base_str[h]))
                    #print bp
                    line_bp += bp
                elif base_str[h].isdigit():
                    bp = int(base_str[h])
                    #print bp
                    line_bp += bp
        i += 1
    return int(line_bp)

total_bp = int(0)

for line in list1:
    line = line.strip()
    line = line.split()
    base_str = line[5]
    line_bp = base_counter(base_str)
    total_bp += line_bp


    
print total_bp
