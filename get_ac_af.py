import myvariant

mv = myvariant.MyVariantInfo()

# read in R output file
filename = "data/allele_ids_to_query.txt"
allele_ids = [line.rstrip('\n') for line in open(filename)]

# remove the header
allele_ids = allele_ids[1:]

# query list
df = mv.querymany(allele_ids, scopes = 'clinvar.allele_id', fields = ['clinvar.allele_id', 'gnomad_exome.af.af', 'gnomad_exome.an.an'], as_dataframe=True)

# filter to only include essential cols
essential_cols = ['gnomad_exome.af.af', 'gnomad_exome.an.an']
df = df[essential_cols]

# remove NaN
df = df.dropna(axis = 0, how='any')

# write df to txt to get it back into R
df.to_csv('data/allele_ids_af_ac.txt')
