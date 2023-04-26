from Bio import Entrez
import ftplib

Entrez.email = "anne.busch@uni-mainz.de"
# handle = Entrez.esearch(db="nucleotide", retmax=10, term="S.cerv[ORGN]")
# handle = Entrez.einfo()#

handle = Entrez.esearch(db="gds", term="cerevisiae, tRNA", retmax=10)
record = Entrez.read(handle)
handle.close()
print(record)


handle = Entrez.efetch(db="gds", id="200198271", retmode="text")
values = handle.readlines()
for entry in values:
    print(entry.strip())
handle.close()

print(values[0])

ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()


# navigate to the directory where the series data is located
ftp.cwd('/geo/series/GSE198nnn/GSE198271/')

# list the files in the directory
filenames = ftp.nlst()
print(filenames)
matrix_filename = [name for name in filenames if name.endswith('matrix')][0]
print(matrix_filename)
ftp.sendcmd("Type I")
size = ftp.size(matrix_filename)
print(size)
# download each file
for matrix_filename in filenames:
    print(matrix_filename)
    # open a local file to write the downloaded data to
    localfile = open(f"{matrix_filename}.txt", 'w')
    # use RETR command to download the file
    # set the FTP transfer mode to ASCII/text
    # ftp.sendcmd('TYPE A')
    file = ftp.retrlines("RETR " + matrix_filename)
    print(file)
    ftp.retrlines(f"RETR {matrix_filename}", localfile.write)

    # close the local file
    localfile.close()

# close the FTP connection
ftp.quit()

