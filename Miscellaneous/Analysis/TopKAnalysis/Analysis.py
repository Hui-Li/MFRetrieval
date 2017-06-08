#-*-coding:utf-8-*-
from __future__ import division

if __name__=="__main__":
    data = "Netflix"
    k = 50
    innerProduct = {}

    for i in xrange(1, k+1):
        innerProduct[i] = 0

    count = 0
    with open(data + "-50.txt") as input:

        currentUserID = None
        currentUserResult = []
        for line in input.readlines():
            line = line.strip()
            if line == "":
                break

            count = count + 1

            paras = line.split(" ")
            if currentUserID == None:
                currentUserID = paras[0]
            elif paras[0]!=currentUserID:
                currentUserResult.sort(reverse=True)
                for i in xrange(1, k+1):
                    innerProduct[i] = innerProduct[i] + currentUserResult[i-1]

                currentUserResult = []
                currentUserID = paras[0]

            currentUserResult.append(float(paras[2]))


    count = count / k
    print "count:" + str(count)

    for i in xrange(1, k+1):
        innerProduct[i] = innerProduct[i]/count

    with open(data + "_topk_innerproduct.txt", "w") as output:
        for i in xrange(1, k + 1):
            output.write(str(i) + "," + str(innerProduct[i]) + "\n")




