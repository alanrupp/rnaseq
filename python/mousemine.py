#!/usr/bin/python2

'''
Get ontology information from Intermine API for a list of genes. Genes should be
in a list format.
'''

def query_intermine(genes):
    genes = ', '.join(genes)
    from intermine.webservice import Service
    service = Service("http://www.mousemine.org/mousemine/service")
    query = service.new_query("OntologyAnnotation")
    query.add_constraint("ontologyTerm", "MPTerm")
    query.add_constraint("subject", "SequenceFeature")
    query.add_view(
    "subject.primaryIdentifier", "subject.symbol",
    "subject.sequenceOntologyTerm.name", "ontologyTerm.identifier",
    "ontologyTerm.name", "evidence.publications.pubMedId",
    "evidence.comments.type", "evidence.comments.description"
    )
    query.add_sort_order("OntologyAnnotation.ontologyTerm.name", "ASC")
    query.add_constraint("subject.organism.taxonId", "=", "10090", code = "A")
    query.add_constraint("subject", "LOOKUP", genes, code = "B")
    query.outerjoin("evidence.comments")
    return query

# convert list of unicode to list of string
def to_string(list_object):
    list_object = [a.encode('ascii', 'ignore') for a in list_object]
    return list_object

# grab relevant fields
def clean_result(query):
    result = {'id': [], 'ontology': []}
    for row in query.rows():
        result['id'].append(row['subject.symbol'])
        result['ontology'].append(row["ontologyTerm.name"])
    for key in result.keys():
        result[key] = to_string(result[key])
    return result
