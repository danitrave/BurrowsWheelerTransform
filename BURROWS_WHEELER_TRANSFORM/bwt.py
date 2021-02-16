##This algorithm is the implementation of the one for a  Burrows-Wheeler transform that allows to convert the genome sequence in a compact string.
##The algorithm must be run on the terminal using python version 3.7 to have a more clear and user-friendly visualization.
##The program asks the user to insert a genome string. Then the algorithm will create the Burrows-Wheeler transform (BWT). In the BWT we see the presence
##of an additional character "$" that is insert by the program, marking the end of the sequence (the character is not insered by the user). Once done that, the program
##will use the BWT to reverse the string and re-obtain the initial string that the user has insert by input. If the messaged visualized is YES, the reverse
##is equal to the initial string. If, otherwise, is NO the reverse is different from the initial string. Then the algorithm will ask the user to insert
##a substring in a way to verify its presence and possibly its number of occurences. If present, the program will report the number of occurencs and its
##range (initial and final index in which the subsequence is found) in the sequence. If, otherwise, the subsequence is not present,
##the program will confirm its absence in the genome sequence.  

matrixRot = []   #all these variable are global since they are used and updated by different functions 
F = ''               #last colum containing the last character of each rotation (so our BWT) used for reverse formation and string matching 
L = ''               #first colum containing the first character of each rotation used for reverse formation and string maching
lastCol = []         #same as the two previous ones but seen as a list of tuples reporting the character and the position 
firstCol = []        #for each kind of character. So if A occurs 3 times, the number associated to each A will be: 0,1,2 or 3
array = tuple()      #array reporting offsets 

def ranker(volatile,i):
    if i == '$':     #if the last character is the special one
        output = '$'     #just report it
        return output
    else:              #if the last character is a string charcter (not the special one)
        output = volatile.count(i) - 1    #count the number of that kind of character in the slice selected -1 since indexes start from 0 and not 1
        return output

def suffix(M):       #this function calculates the offset of each substring of the input string reported as first character from the special charater "$" from each cyclic rotation
    global array
    for t in M:      #select the rotation from the matrix 
        i = 1
        suffix = len(t)       #take the length of each rotation
        while True:
            if t[i-1] == '$':       #if the charatcer is the special one
                suffix = suffix - i     #calculate offset subtracting the current positon to the total length of the cyclic rotation length 
                array += (suffix,)         #add the result to the array 
                break                  #stop execution when "$" is found 
            else:
                i += 1       #keep on elongating the offset as the special character is not found 
    return array         #return the array containing all offsets

def cyclicRotations(T):   #this funciton use the sequence to build the cyclic rotaions
    global matrixRot
    TT = T * 2
    matrixRot = [TT[i:i+len(T)] for i in range(0,len(T))]  #built cyclic rotations
    matrixRot = sorted(matrixRot)        #sort the rotations in lexycographic order
    array = suffix(matrixRot)      #use the matrix reporting all cyclic rotations as strings to build an array reporting the offset of each rotation
    return matrixRot        #Burrows-Wheeler matrix (BWM)

def BWT(T):   #this function computes the BWT
    global F
    global L
    T = T + '$'    #add special character at the end of the string
    rotations = cyclicRotations(T)    #build cyclic rotations
    for e in rotations:     #of each cyclic rotation take the first character and the last one to get
        L += e[-1]        #the last column of the matrix, so our BWT
        F += e[0]            #the first column
    return L

def firstColRank(F,LL):  #this function ranks the characters of the first column of the BWM
    global firstCol        #use the first column. We can do in this way since characters of a kind of the first column
    for y in F:             #have the same rank-order of the character of a kind in the last column --> LM-mapping rule
        j = 0
        while True:
            currentCell = tuple()      #tuple that will collect the character and the rank
            if LL[j][0] == y:            #take the first element of the last column-ranked and the first one of the first column not ranked
                currentCell = (y,) + (LL[j][1],)   #associate the same rank of the first character of a kind of the last column to the first character of the first column of kind
                firstCol += (currentCell,)    #build up the new character-rank first column
                LL.remove(LL[j])         #remove the character of the last column that we have seleced to avoid ambiguos and repetitive ranking of the first column
                break              #stop while loop
            else:
                j += 1      #go to the next character 
    return firstCol
                    
def T_ranking():       #this function builts the ranks of the first and last column of the BWM
    global lastCol
    ranking_L = []      #collector for the ranks of the last characters 
    i = 0
    for v in matrixRot:
        slicer = v.find('$')   #find the index of the special character at each rotation
        volatile = v[slicer:]      #take the section fo each rotation starting from the special charater up to the end of the string
        currentRank = ranker(volatile,L[i])    #rank the last character of each cyclic rotation in the matrix
        ranking_L += [currentRank,]     #add the current 
        i += 1
    i = 0
    while i != len(L):    #build the list of tuples reporting the character and the rank of that character 
        currentCell = (L[i],) + (ranking_L[i],)     #create the current character-rank tuple
        lastCol += (currentCell,)                #add it to the the new lst column
        i += 1
    LL = lastCol[:]                  #make a copy to avoid unwanted alliasing effects of the last column-rank
    firstCol = firstColRank(F,LL)       #build rank of the characters of the first column using the last column (character-rank) as reference
    return firstCol,lastCol       #report the first and the last columns each seen as a list of charcter-rank tuples 

def reverseBWT(T):         #this function reverse the BWT to get the original strinf
    colums = T_ranking()       #use first as a sort of guide line and last column to get the reverse 
    reverse = ''
    guard = colums[1][0]          #the starting point is the first character of the last column
    while guard != ('$','$'):        #until you do not have the special character
        reverse += guard[0]             #add current character of the last character 
        index = colums[0].index(guard)    #use the first column as an array to get the index of the column having the same character of the last column  
        guard = colums[1][index]      #move to the character of the last column at the index that we just found before. 

    return reverse[::-1]          #get reverse 

def stringMatch(P):       #This functions takes a string as parameter and returns the number of its occurences in the string and the range of the match
    global firstCol
    global lastCol
    global array

    i = len(P)-1       #start from the last character of the string i want to match 
    start = []
    suffix = P[i]                    #the last character of the string i want to match is my new suffix 
    for x in range(len(firstCol)):       #consider all elements in the ranked last column 
        vol = tuple()
        if firstCol[x][0] == suffix:        #if the character of the last column-ranked is equal to our suffix 
            vol = (firstCol[x],) + (x,)        #collect it. It will be one of the possible starting points from which i start to match the remainig characters of my string 
            start += (vol,)
    
    if len(P) == 1:            #in the case in which the substring is rapresented by one character only 
        indexForRange = []           #the occurrences are present in the initials point extrapolated from the first column-ranked  
        count = 0                 
        if len(start) == 0:             #of course i evaluate the occurences of the unique character only if it is present the the string
            return count,indexForRange
        else:
            for char in start:           #if it is present at least one time in the string
                count += 1                        #i collect the number of occurences of it 
                index = firstCol.index(char[0])         #and i collect the occurence position the the BWM 
                indexForRange += [index,]                  #to then use it to extrapolate the range of the occurences
    
    else:                        #if, otherwise, the substring is made of more then one character
        count = 0
        indexForRange = []           #initialize a list that will contain the ranges of the occurences of the substring
        for e in start:                #start form the first possible element that may lead to a complete match of th substring
            i = len(P)-1                 #set as suffix the slice of the substring: as long as a match is found, the slice is increased un to get the whole sequence when its value is 0
            index = e[1]                   #set the value of the first current colomn-rank
            while True:
                i = i - 1                    #move to the next character of the substring. The current suffix length will be increased
                guard = lastCol[index]                #select the character of the last column-ranked localized in the same index as the character of the first column-rank
                if guard[0] == P[i]:                     #if the character of the tuple of last column-ranked is equal to the current character of the suffix
                    index = firstCol.index(guard)             #go on to the new character
                    if i == 0:                                   #when the suffix reaches the max lenth, means that you found a match
                        indexForRange += [index,]           #select the index of the BWM where the last character of the substring has been matched
                        count += 1                  #keep on track of the number of matches you found
                        break
                else:
                    break             #if no match is found, go to next starting point 
            
    rangeSub = []                  #list that will be filled with lists each reporting the start-end index of each match per substring found 
    for r in indexForRange:
        rank = []                    #range of the current match. The range is evaluated as python indexes. So the last charater is index -1
        loc = array[r]                   #use the array to report the offset of a string in a given position 
        rank = [loc,(loc+len(P)-1)]         #select start and end of the occurence 
        rangeSub += [rank,]           #add current range 

    return count,sorted(rangeSub)     #report number of occurences and for each, the range in the string 

T = input('Insert the string to get its BWT:')
T = T.upper()
BWT = BWT(T)
print('The BWT obtained is:',BWT)
reverse = reverseBWT(T)
print('Using the BWT, the reverse string obtained is:',reverse)
if reverse == T:
    result = True
else:
    result = False
if result == True:
    print('Is the reverse string that we got equal to the original string? YES')
else:
    print('Is the reverse string that we got equal to the original string? NO')

subString = input('Insert the substring for which you want to find the occurence/s:')
subString = subString.upper()
match = stringMatch(subString)
if match[0] != 0:
    print('The substring occurs in the string',match[0],'times at index/s',match[1])
else:
    print('The substring insered is not present in the string')

    
    
