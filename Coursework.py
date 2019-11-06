def dynprog(letters,scoring,seq1,seq2):
    table = [["-","-"],["-"]]
    for i in range(len(seq1)):
        table[0].append(seq1[i])
    for i in range(len(seq2)):
        table.append([seq2[i]])
    for i in range(1,len(table)):
        for j in range(len(seq1)+1):
            table[i].append(0)
    directions = [["-","-"],["-"]]
    for i in range(len(seq1)):
        directions[0].append(seq1[i])
    for i in range(len(seq2)):
        directions.append([seq2[i]])
    for i in range(1,len(directions)):
        for j in range(len(seq1)+1):
            directions[i].append([])
    directions[1][1]=["f"]
    for i in range(2,len(seq2)+2):
        for j in range(len(letters)):
            if letters[j] == seq2[i-2]:
                table[i][1] = table[i-1][1]+max(0,int(scoring[j][-1]))
        directions[i][1] = ["u"]
    for i in range(2,len(seq1)+2):
        for j in range(len(letters)):
            if letters[j] == seq1[i-2]:
                table[1][i] = table[1][i-1]+max(0,int(scoring[j][-1]))
        directions[1][i] = ["l"]
    for i in range(2,len(seq1)+2):
        for j in range(2,len(seq2)+2):
            letter1 = table[0][i]
            letter2 = table[j][0]
            score = 0
            indel1 = 0 
            indel2 = 0
            for k in range(len(letters)):
                for l in range(len(letters)):
                    if letters[k].upper() == letter1.upper() and letters[l].upper()==letter2.upper():
                        score = int(scoring[k][l])
                        indel1 = int(scoring[k][-1])
                        indel2 = int(scoring[-1][l])
            table[j][i] = max(0,table[j-1][i-1]+score,table[j-1][i]+indel2,table[j][i-1]+indel1)
            if table[j][i] == table[j-1][i-1]+score:
                directions[j][i]="d"
            elif table[j][i] == table[j-1][i]+indel2:
                directions[j][i]="u"
            elif table[j][i] == table[j][i-1]+indel1:
                directions[j][i]="l"
    max_num = 0
    x = 1
    y = 1
    for j in range(1,len(table)):
        for i in range(1,len(table[0])):
            if table[j][i] >= max_num:
                max_num = table[j][i]
                x = i
                y = j
    new_seq1,new_seq2 = local_alignment(table,directions,x,y)
    new_seq1 = reverse(new_seq1)
    new_seq2 = reverse(new_seq2)
    indices1,indices2 = matches(new_seq1,new_seq2)
    return [max_num,indices1,indices2]

def dynproglin(letters,scoring,seq1,seq2):
    new_seq1,new_seq2,score = linear(letters,scoring,seq1,seq2)
    new_seq1 = reverse(new_seq1)
    new_seq2 = reverse(new_seq2)
    indices1,indices2 = matches(new_seq1,new_seq2)
    return [score,indices1,indices2]

def heuralign(letters,scoring,seq1,seq2):
    if min(len(seq1),len(seq2))>10:
        ktup = 3
    elif min(len(seq1),len(seq2))>5:
        ktup =2
    else:
        ktup = 1
    index = {}
    if len(seq1)>=len(seq2):
        long_seq = seq1
        short_seq = seq2
    else:
        long_seq = seq2
        short_seq = seq1
    for i in range(len(long_seq)-ktup+1):
        if long_seq[i:i+ktup] in index:
            index[long_seq[i:i+ktup]].append(i)
        else:
            index[long_seq[i:i+ktup]]=[i]
    pairs = []
    dists = []
    for i in range(len(short_seq)-ktup+1):
        if short_seq[i:i+ktup] in index:
            if len(index[short_seq[i:i+ktup]])>0:
                for j in range(len(index[short_seq[i:i+ktup]])):
                    pairs.append([index[short_seq[i:i+ktup]][j],i])
                    dists.append(pairs[-1][0]-i)
    potentials = []
    scores = []
    locations = []
    for i in range(len(dists)):
        for j in range(i,len(dists)):
            if dists[i]==dists[j]:
                alignment1 =long_seq[pairs[i][0]:pairs[i][0]+ktup]
                alignment2 = short_seq[pairs[i][1]:pairs[i][1]+ktup]
                score = calc_score(letters,scoring,alignment1,alignment2)
                new_alignment1 = long_seq[pairs[i][0]:pairs[j][0]+ktup]
                new_alignment2 = short_seq[pairs[i][1]:pairs[j][1]+ktup]
                if len(new_alignment1) == len(new_alignment2):
                    new_score = calc_score(letters,scoring,new_alignment1,new_alignment2)
                    if new_score>=score:
                        potentials.append([new_alignment1,new_alignment2])
                        scores.append(new_score)
                        locations.append([pairs[i][0],pairs[i][1],len(new_alignment1)])
                    else:
                        potentials.append([alignment1,alignment2])
                        scores.append(score)
                        locations.append([pairs[i][0],pairs[i][1],len(alignment1)])
    if len(scores)==0:
        return [0,[],[]]
    max_score = max(scores)
    best = []
    best_locations = []
    for i in range(len(potentials)):
        if scores[i]==max_score or scores[i] == max_score-1:
            best.append(potentials[i])
            best_locations.append(locations[i])
    points = []
    for i in range(len(best_locations)):
        for j in range(best_locations[i][2]):
            if long_seq==seq1:
                points.append([best_locations[i][0]+j,best_locations[i][1]+j,i])
            else:
                points.append([best_locations[i][1]+j,best_locations[i][0]+j,i])
    new_sequences = []
    for i in range(len(points)):
        for k in range(len(points)):
            if [points[i][0]+1,points[i][1],k] in points:
                for j in range(len(points)):
                    if points[i][2] == points[j][2]:
                        start = points[j]
                        break
                for j in range(len(points)):
                    if points[j][0] == points[i][0]+1 and points[j][1] == points[i][1]:
                        k = points[j][2]
                        break
                for j in range(len(points)):
                    if points[j][2] == k:
                        end = points[j]
                new_sequences.append([seq1[start[0]:end[0]],seq2[start[1]:points[i][1]]+"-"+seq2[points[i][1]:end[1]]])
            if [points[i][0],points[i][1]+1,k] in points:
                for j in range(len(points)):
                    if points[i][2] == points[j][2]:
                        start = points[j]
                        break
                for j in range(len(points)):
                    if points[j][0] == points[i][0] and points[j][1] == points[i][1]+1:
                        k = points[j][2]
                        break
                for j in range(len(points)):
                    if points[j][2] == k:
                        end = points[j]
                new_sequences.append([seq1[start[0]:points[i][0]]+"-"+seq1[points[i][0]:end[0]],seq2[start[1]:end[1]]])
            if [points[i][0]+2,points[i][1],k] in points:
                for j in range(len(points)):
                    if points[i][2] == points[j][2]:
                        start = points[j]
                        break
                for j in range(len(points)):
                    if points[j][0] == points[i][0]+2 and points[j][1] == points[i][1]:
                        k = points[j][2]
                        break
                for j in range(len(points)):
                    if points[j][2] == k:
                        end = points[j]
                new_sequences.append([seq1[start[0]:end[0]],seq2[start[1]:points[i][1]]+"--"+seq2[points[i][1]:end[1]]])
            if [points[i][0],points[i][1]+2,k] in points:
                for j in range(len(points)):
                    if points[i][2] == points[j][2]:
                        start = points[j]
                        break
                for j in range(len(points)):
                    if points[j][0] == points[i][0] and points[j][1] == points[i][1]+2:
                        k = points[j][2]
                        break
                for j in range(len(points)):
                    if points[j][2] == k:
                        end = points[j]
                new_sequences.append([seq1[start[0]:points[i][0]]+"--"+seq1[points[i][0]:end[0]],seq2[start[1]:end[1]]])
    for i in range(len(new_sequences)):
        best.append(new_sequences[i])
    max_score =0
    best_index = 0
    for i in range(len(best)):
        score = calc_score(letters,scoring,best[i][0],best[i][1])
        if score>=max_score:
            best_index = i
            max_score = score
    new_seq1 = best[best_index][0]
    new_seq2 = best[best_index][1]
    indices1,indices2 = matches(new_seq1,new_seq2)
    return [max_score,indices1,indices2]

def linear(letters,scoring,seq1,seq2):
    new_seq1 = ""
    new_seq2 = ""
    table = [[0],[0]]
    max_score = 0
    x=0
    y=0
    direction = "0"
    letter1 = "-"
    letter2 = "-"
    for i in range(1,len(seq1)+1):
        for j in range(len(letters)):
            if letters[j] == seq1[i-1]:
                table[0].append(max(table[0][i-1]+int(scoring[j][-1]),0))
        table[1].append(0)
    if len(seq2)==0:
        return "","",0
    for i in range(len(seq2)):
        for j in range(len(letters)):
            if letters[j]==seq2[i]:
                table[1][0]=table[0][0]+max(0,int(scoring[j][-1]))
        for j in range(1,len(seq1)+1):
            score = 0
            indel1 = 0 
            indel2 = 0
            for k in range(len(letters)):
                for l in range(len(letters)):
                    if letters[k].upper() == seq1[j-1].upper() and letters[l].upper()==seq2[i].upper():
                        score = int(scoring[k][l])
                        indel1 = int(scoring[k][-1])
                        indel2 = int(scoring[-1][l])
            diagonal = table[0][j-1] +score
            up = table[0][j]+indel2
            left = table[1][j-1]+indel1
            table[1][j]=(max(diagonal,up,left,0))
            if table[1][j]>max_score:
                max_score = table[1][j]
                x=j
                y=i
                if max_score == table[0][j-1]+score:
                    direction = "d"
                elif max_score == table[0][j]+indel2:
                    direction="u" 
                elif max_score == table[1][j-1]+indel1:
                    direction = "l"
                letter1 = seq1[j-1]
                letter2 = seq2[i]
        for m in range(len(table[1])):
            table[0][m] = table[1][m]
    if direction == "d":
        new_seq1 = new_seq1+letter1
        new_seq2 = new_seq2 +letter2
        x-=1
        y-=1
    elif direction == "u":
        new_seq1 = new_seq1+"-"
        new_seq2 = new_seq2 +letter2
        y-=1
    elif direction == "l":
        new_seq1 = new_seq1+letter1
        new_seq2 = new_seq2+"-"
        x-=1
    while True:
        table=[[0],[0]]
        for i in range(1,len(seq1)+1):
            for j in range(len(letters)):
                if letters[j] == seq1[i-1]:
                    table[0].append(max(table[0][i-1]+int(scoring[j][-1]),0))
            table[1].append(0)
        if y+1>0:
            for i in range(y+1):
                for j in range(len(letters)):
                    if letters[j]==seq2[i]:
                        table[1][0]=table[0][0]+max(0,int(scoring[j][-1]))
                if i==y and x==0 and table[1][0] ==0:
                    return new_seq1,new_seq2,max_score
                for j in range(1,len(seq1)+1):
                    score = 0
                    indel1 = 0 
                    indel2 = 0
                    for k in range(len(letters)):
                        for l in range(len(letters)):
                            if letters[k].upper() == seq1[j-1].upper() and letters[l].upper()==seq2[i].upper():
                                score = int(scoring[k][l])
                                indel1 = int(scoring[k][-1])
                                indel2 = int(scoring[-1][l])
                    diagonal = table[0][j-1] +score
                    up = table[0][j]+indel2
                    left = table[1][j-1]+indel1
                    table[1][j]=(max(diagonal,up,left,0))
                    if i==y and j==x:
                        if table[1][j] == 0:
                            return new_seq1,new_seq2,max_score
                        elif table[1][j] == table[0][j-1]+score:
                            direction = "d"
                        elif table[1][j] == table[0][j]+indel2:
                            direction="u"
                        elif table[1][j] == table[1][j-1]+indel1:
                            direction = "l"
                        letter1 = seq1[j-1]
                        letter2 = seq2[i]
                for m in range(len(table[1])):
                    table[0][m] = table[1][m]
        else:
            for i in range(len(seq2)):
                for j in range(len(letters)):
                    if letters[j]==seq2[i]:
                        table[1][0]=table[0][0]+max(0,int(scoring[j][-1]))
            if table[0][x] == 0:
                return new_seq1,new_seq2,max_score
            else:
                direction = "l"
                new_seq1 = new_seq1+letter1
                new_seq2 = new_seq2+"-"
                x-=1
        if direction == "d":
            new_seq1 = new_seq1+letter1
            new_seq2 = new_seq2 +letter2
            x-=1
            y-=1
        elif direction == "u":
            new_seq1 = new_seq1+"-"
            new_seq2 = new_seq2 +letter2
            y-=1
        elif direction == "l":
            new_seq1 = new_seq1+letter1
            new_seq2 = new_seq2+"-"
            x-=1

def calc_score(letters,scoring,seq1,seq2):
    score = 0
    indel1 = 0 
    indel2 = 0
    for i in range(len(seq1)):
        found = False
        letter1 = seq1[i]
        letter2 = seq2[i]
        for k in range(len(letters)):
            for l in range(len(letters)):
                if letters[k].upper() == letter1.upper() and letters[l].upper()==letter2.upper():
                    score += int(scoring[k][l])
                    found = True
        if found == False:
            for k in range(len(letters)):
                if letters[k].upper() ==  letter1.upper():
                    score += int(scoring[k][-1])
                    found = True
                    break
        if found == False:          
            for l in range(len(letters)):
                if letters[l].upper() == letter2.upper():
                    score += int(scoring[-1][l])
                    break
    return score

def matches(new_seq1,new_seq2):
    indices1 = []
    indices2 = []
    for i in range(len(new_seq1)):
        if new_seq1[i]!="-" and new_seq2[i]!="-":
            counter = 0
            for j in range(i):
                if new_seq1[j]=="-":
                    counter+=1
            indices1.append(i-counter)
            counter = 0
            for j in range(i):
                if new_seq2[j]=="-":
                    counter+=1
            indices2.append(i-counter)
    return indices1,indices2

def reverse(sequence):
    return sequence[::-1]

def local_alignment(table,directions,x,y):
     new_seq1 = ""
     new_seq2 = ""
     while table[y][x]!=0:
         if directions[y][x] == "d":
             new_seq1 = new_seq1+table[0][x]
             new_seq2 = new_seq2+table[y][0]
             x-=1
             y-=1
         elif directions[y][x] == "u":
             new_seq1=new_seq1+"-"
             new_seq2=new_seq2+table[y][0]
             y-=1
         else:
             new_seq1=new_seq1+table[0][x]
             new_seq2=new_seq2+"-"
             x-=1
     return new_seq1,new_seq2
