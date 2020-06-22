# Finding the sets of attributed that are accessed for the most part, by distinct sets of applications (Queries)
# We look for dividing points along the diagonal such that
# 1. Total accesses to only onr fragment is maxmized, while
# 2. Total accesses to more tham one fragments are minimized
import BEA
import math

QA_matrix = BEA.query_access_matrix  # Query access matrix
sum_attr_access = BEA.sum_attr_access  # Sum of accesses by the query applications
order_CA = BEA.CA[0, :]
order_CA = order_CA.astype(int)
no_of_queries = BEA.no_of_queries


# Taking the sums of accesses of queries
def sum_access(arr):
    sum = 0
    for items in arr:
        sum = sum + sum_attr_access[items - 1]
    return sum


# Find the partitioning point such that cost Z is maximized
def partition():
    z = []
    fragments = []
    for split_point in range(1, len(order_CA)):
        frag1 = order_CA[0:split_point]
        frag2 = order_CA[split_point:len(order_CA)]
        fragments.append([frag1, frag2])
        print("Fragments =", frag1, frag2)

        TA = []  # Cluster for attributes in frag 1
        TB = []  # Cluster for attributes in frag 2
        for i in range(no_of_queries):
            use_frag1 = 0
            for items in frag1:
                # if queries accesses any of the attributes in fragment 1
                # then TA of that query will be 1
                if QA_matrix[i, items - 1] == 1:
                    use_frag1 = 1
                    break
            TA.append(use_frag1)

            use_frag2 = 0
            for items in frag2:
                # if queries accesses any of the attributes in fragment 2
                # then TB of that query will be 1
                if QA_matrix[i, items - 1] == 1:
                    use_frag2 = 1
                    break
            TB.append(use_frag2)
        print("TA =", TA)
        print("TB =", TB)

        # Queries are classified only into three sets
        TQ = []  # Applications that only use attributes in TA
        BQ = []  # Applications that only use attributes in BA
        OQ = []  # Applications that only use attributes in both TA and BA
        for i in range(no_of_queries):
            if TA[i] == 0 and TB[i] == 1:
                BQ.append(i + 1)
            elif TA[i] == 1 and TB[i] == 0:
                TQ.append(i + 1)
            else:
                OQ.append(i + 1)
        print("TQ =", TQ)
        print("BQ =", BQ)
        print("OQ =", OQ)

        CTQ = sum_access(TQ)  # Total number of accesses to attributes by TQ
        CBQ = sum_access(BQ)  # Total number of accesses to attributes by BQ
        COQ = sum_access(OQ)  # Total number of accesses to attributes by OQ

        # Find the partitioning point such that cost Z is maximized
        z.append(CTQ * CBQ - math.pow(COQ, 2))
        print("z =", z)

    if max(z) < 0:
        print("Vertical Fragmentation not possible.")
    else:
        print("Best partitioning point = ", fragments[z.index(max(z))][0], fragments[z.index(max(z))][1])


partition()
