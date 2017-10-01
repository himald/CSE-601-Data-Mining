import itertools
import copy


def get_gene_matrix(file_name):
    # store the gene expression file
    file_handle = open(file_name + ".txt", "r")
    lines = file_handle.readlines()

    # initialize gene matrix
    gene_matrix = []
    gene_matrix_flattened_list = []

    # read each line in file
    for line in lines:
        # split line and store as expression list
        expression_list = line.strip().split("\t")

        # append G_ to all the expressions and exclude the disease names
        gene_expression_list = []
        j = 1
        for expression in expression_list:
            if j == len(expression_list):
                gene_expression_list.append(expression)
            else:
                gene_expression_list.append("G" + str(j) + "_" + expression)
            j += 1

        # populate the gene expression list in gene matrix and flattened gene matrix list
        gene_matrix.append(gene_expression_list)
        gene_matrix_flattened_list += gene_expression_list
    return gene_matrix, gene_matrix_flattened_list


def get_frequent_itemsets(gene_matrix, gene_matrix_flattened_list, support_threshold_percentage):
    # initialize support threshold
    support_threshold = support_threshold_percentage / 100

    # store the unique elements in flattened gene matrix list
    gene_matrix_unique_list = list(set(gene_matrix_flattened_list))

    # declare 2 matrices, one for master dictionary and one for each iteration
    master_gene_dict = dict()
    current_iter_dict = dict()

    # copy the gene_matrix_unique_list for iteration
    itemset_candidate_list = gene_matrix_unique_list

    total_entries = len(gene_matrix)

    # obtain combinations of length 1 to size of gene_matrix_unique_list until none of the combinations exceed minimal support
    for itemset_length in range(1, len(gene_matrix_unique_list) + 1):
        # initialize the master_gene_dict and current_iter_dict
        initialize_dictionaries(itemset_length, itemset_candidate_list, current_iter_dict, master_gene_dict)

        # populate the support for the combinations
        populate_itemset_support(gene_matrix, master_gene_dict, current_iter_dict, total_entries)

        # reinitialize current iteration dictionary and itemset_candidate_list
        current_iter_dict = dict()
        itemset_candidate_list = []

        # populate itemset_candidate_list with itemsets which have support greater than the minimal support threshold
        for itemset_string in master_gene_dict:
            if (master_gene_dict[itemset_string] >= support_threshold and str(itemset_string).count(
                    ",") == itemset_length - 1):
                itemset_candidate_list.append(itemset_string)

        # break the loop if itemset_candidate_list is empty
        if len(itemset_candidate_list) == 0:
            break

    frequent_itemsets_dict = dict()
    for itemset in master_gene_dict:
        itemset_support = master_gene_dict[itemset]
        if itemset_support >= support_threshold:
            itemset_length = len(itemset.split(","))
            if itemset_length not in frequent_itemsets_dict:
                frequent_itemsets_dict[itemset_length] = []
            frequent_itemsets_dict[itemset_length].append(itemset)

    return frequent_itemsets_dict


def initialize_dictionaries(itemset_length, itemset_candidate_list, current_iter_dict, master_gene_dict):
    # for 1 combination itemsets, initialize the master_gene_dict and current_iter_dict elements to 0
    if itemset_length == 1:
        for itemset in itemset_candidate_list:
            master_gene_dict[itemset] = 0
            current_iter_dict[itemset] = 0

    # for 2 or more combination itemsets, form new itemsets by taking union of existing itemsets
    else:
        for i in range(0, len(itemset_candidate_list)):
            itemset_A = itemset_candidate_list[i]
            itemset_A_list = itemset_A.split(",")
            itesmset_A_list_set = set(itemset_A_list)
            for j in range(i + 1, len(itemset_candidate_list)):
                itemset_B = itemset_candidate_list[j]
                itemset_B_list = itemset_B.split(",")
                itemset_B_list_set = set(itemset_B_list)

                # check if both itemsets are same or if they don't have enough common elements
                intersection_length = len(itesmset_A_list_set & itemset_B_list_set)
                if (intersection_length == itemset_length - 1 or intersection_length < itemset_length - 2):
                    continue
                new_itemset_list = list(itesmset_A_list_set | itemset_B_list_set)

                # sort the new itemset to prevent a new key being created for the same permutation
                new_itemset_list.sort()

                # convert the new itemset list to string
                new_itemset_string = ",".join(new_itemset_list)

                # add the new itemset string to both dictionaries
                master_gene_dict[new_itemset_string] = 0
                current_iter_dict[new_itemset_string] = 0


def populate_itemset_support(gene_matrix, gene_dict, current_dict, total_entries):
    for gene_list in gene_matrix:
        gene_list_set = set(gene_list)
        for itemset in current_dict:
            itemset_list = itemset.split(",")

            # if combination is a subset of genes in the row
            if (set(itemset_list) < gene_list_set):
                gene_dict[itemset] += 1

    for item_set in current_dict:
        gene_dict[item_set] = gene_dict[item_set] / total_entries
    return gene_dict


def generate_rules(rules_set, frequent_itemsets_dict, confidence_threshold_percentage, gene_matrix):
    frequent_itemsets_keylist = list(frequent_itemsets_dict.keys())
    frequent_itemsets_keylist.sort(reverse=True)

    for itemset_length in frequent_itemsets_keylist:
        # exclude itemsets of length 1 for rule generation
        if itemset_length == 1:
            break
        else:
            itemset_list = frequent_itemsets_dict[itemset_length]
            for itemset_string in itemset_list:
                itemset_list_set = set(itemset_string.split(","))
                if itemset_list_set not in rules_set:
                    body_max_length = itemset_length - 1
                    numerator = 0
                    for gene_list in gene_matrix:
                        gene_list_set = set(gene_list)
                        if itemset_list_set < gene_list_set:
                            numerator += 1
                    if numerator > 0:
                        itemset_string_list = itemset_string.split(",")
                        get_itemset_rules(itemset_string_list, body_max_length, numerator,
                                          confidence_threshold_percentage, gene_matrix, rules_set)

def get_itemset_rules(itemset_list, body_max_length, numerator, confidence_threshold_percentage, gene_matrix, rules_set):
    itemset_list_set = set(itemset_list)

    for body_length in range(body_max_length, 0, -1):
        for body_list in list(itertools.combinations(itemset_list, body_length)):
            body_set = set(body_list)
            head_list = list(body_set ^ itemset_list_set)
            rule = ",".join(body_list) + "->" + ",".join(head_list)
            if rule not in rules_set:
                denominator = 0
                for gene_list in gene_matrix:
                    gene_list_set = set(gene_list)
                    if body_set < gene_list_set:
                        denominator += 1
                confidence_percentage = (numerator / denominator) * 100
                if (confidence_percentage >= confidence_threshold_percentage):
                    rules_set.add(rule)
    return rules_set


def get_template1_input():
    template1_query = input('Enter template-1 query (RULE|BODY|HEAD;ANY|NUMBER|NONE;ITEM1,ITEM2,...): ').split(';')
    template1_query_part = template1_query[0]
    template1_constraint = template1_query[1]
    template1_genelist = template1_query[2].split(',')
    return template1_query_part, template1_constraint, template1_genelist


def get_template1_rules(rules_set, filtered_rules_set, template3_operator):
    template1_query_part, template1_constraint, template1_genelist = get_template1_input()
    template1_query_part = template1_query_part.upper()
    template1_constraint = template1_constraint.upper()
    for rule in rules_set:
        superset = get_superset(rule, template1_query_part)
        if template1_constraint == 'ANY' and len(set(template1_genelist) & superset) > 0:
            filtered_rules_set.add(rule)
        elif template1_constraint == 'NONE' and len(set(template1_genelist) & superset) == 0:
            filtered_rules_set.add(rule)
        elif template1_constraint == '1' and len(set(template1_genelist) & superset) == 1:
            filtered_rules_set.add(rule)
        elif template3_operator == 'AND' and rule in filtered_rules_set:
            filtered_rules_set.remove(rule)


def get_template2_input():
    template2_query = input('Enter template-2 query (RULE|BODY|HEAD;SIZE): ').split(';')
    template2_query_part = template2_query[0]
    template2_size = int(template2_query[1])
    return template2_query_part, template2_size


def get_template2_rules(rules_set, filtered_rules_set, template3_operator):
    template2_query_part, template2_size = get_template2_input()
    template2_query_part = template2_query_part.upper()
    for rule in rules_set:
        superset = get_superset(rule, template2_query_part)
        if len(superset) >= template2_size:
            filtered_rules_set.add(rule)
        elif template3_operator == 'AND' and rule in filtered_rules_set:
            filtered_rules_set.remove(rule)


def get_template3_input():
    template3 = input('Enter template combination: ')
    template3 = template3.upper()
    if 'OR' in template3:
        template3_first = template3.split('OR')[0]
        template3_second = template3.split('OR')[1]
        template3_operator = 'OR'
    elif 'AND' in template3:
        template3_first = template3.split('AND')[0]
        template3_second = template3.split('AND')[1]
        template3_operator = 'AND'
    
    return template3_first, template3_second, template3_operator


def get_superset(rule, template_query_part):
    rule_list = rule.split('->')
    body_list = rule_list[0].split(',')
    head_list = rule_list[1].split(',')
    if template_query_part == 'RULE':
        superset = set(body_list) | set(head_list)
    elif template_query_part == 'BODY':
        superset = set(body_list)
    elif template_query_part == 'HEAD':
        superset = set(head_list)
    return superset


support_threshold_percentage = int(input('Enter support percentage: '))
confidence_threshold_percentage = int(input('Enter confidence percentage: '))
print("Support is set to be " + str(support_threshold_percentage) + "%")


# fetch the gene matrix and its flattened list from the file
gene_matrix, gene_matrix_flattened_list = get_gene_matrix("associationruletestdata")

# get the frequent itemsets
frequent_itemsets_dict = get_frequent_itemsets(gene_matrix, gene_matrix_flattened_list, support_threshold_percentage)

# print the number of frequent itemsets
total_itemsets = 0
for itemset_length in frequent_itemsets_dict:
    num_of_itemsets = len(frequent_itemsets_dict[itemset_length])
    print("Number of length-" + str(itemset_length) + " frequent itemsets: " + str(num_of_itemsets))
    total_itemsets += num_of_itemsets

print("Number of all lengths frequent itemsets: " + str(total_itemsets))


rules_set = set()
generate_rules(rules_set, frequent_itemsets_dict, confidence_threshold_percentage, gene_matrix)
print("Number of rules for support " +str(support_threshold_percentage)+ "% and confidence " +str(confidence_threshold_percentage)+ "% : "+str(len(rules_set)))
template = input('Enter template type: ')
filtered_rules_set = set()
template3_operator = ''
if template == '1':
    get_template1_rules(rules_set, filtered_rules_set, template3_operator)

elif template == '2':
    get_template2_rules(rules_set, filtered_rules_set, template3_operator)

elif template == '3':
    template3_first, template3_second, template3_operator = get_template3_input()
    if template3_first == '1':
        get_template1_rules(rules_set, filtered_rules_set, template3_operator)
    elif template3_first == '2':
        get_template2_rules(rules_set, filtered_rules_set, template3_operator)
    if template3_operator == 'AND':
        updated_rules_set = copy.deepcopy(filtered_rules_set)
    elif template3_operator == 'OR':
        updated_rules_set = copy.deepcopy(rules_set)
    if template3_second == '1':
        get_template1_rules(updated_rules_set, filtered_rules_set, template3_operator)
    elif template3_second == '2':
        get_template2_rules(updated_rules_set, filtered_rules_set, template3_operator)




print(filtered_rules_set)
print("Count of filtered rules: " + str(len(filtered_rules_set)))