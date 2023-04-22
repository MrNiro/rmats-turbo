import sys
import os
import subprocess


def prepare_data(filepath):
    for _, _, filenames in os.walk(filepath):
        for filename in filenames:
            if "MATS.JCEC.txt" in filename:
                all_genes = {}
                f = open(filepath + filename)
                new_filename = filename.replace(".", "_", 2)
                new_filename = new_filename.replace("txt", "csv")
                g = open("./results/" + new_filename, "w")
                g.write("gene_id,gene_name,p_value\n")
                for line in f.readlines()[1:]:
                    info = line.strip().split("\t")
                    gene_name = info[1].replace("\"", "")
                    gene_id = gene_name
                    p_value = info[-5]
                    if gene_name not in all_genes:
                        all_genes[gene_name] = p_value
                        # print(gene_name, p_value)
                        g.write(gene_id + "," + gene_name + "," + p_value + "\n")
                f.close()
                g.close()


def pathway_enrichment():
    for _, _, filenames in os.walk("./results"):
        for filename in filenames:
            if "csv" in filename:
                cmd = "Rscript pathway_enrich.R ./results/" + filename
                print(cmd)
                status, output = subprocess.getstatusoutput(cmd)
                if status != 0:
                    print("\n\nError when enrichment for:", filename)
                    print("Error Message:")
                    print(output)
                    print("\n\n")


if __name__ == '__main__':
    file_path = "D:/Lin/非我有/Harvard Medical School/rMATs/results-1026-IFN/"
    if file_path == "":
        if len(sys.argv) > 1:
            file_path = sys.argv[1]
        else:
            raise Exception("Please provide a valid path to rMATS results!")
    prepare_data(file_path)

    pathway_enrichment()
