using DataFrames
using CSV
using VariantCallFormat
# import Pkg; Pkg.add("DataFrames")

file = "/home/annebusch/Documents/anne/PhD/bedMethyl/modbam2bed-master/test_data/6mA.csv"
content = CSV.read(file, DataFrame, delim="\t", header=false)


reader = VCF.Reader(open("/home/annebusch/Documents/anne/PhD/pinf_sc50.vcf", "r"))
records = read(reader)
print(records)
close(reader)


record = VCF.Record("20\t145430\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
VCF.pos(record)

header = reader.header
print(header)

# iterate over header using metainfo indices
print(header.metainfo[6])
for i in header.metainfo
    println(i)
end


# make new metainfo entry like this (requires dummy metainfo to overwrite)
metainfo = VCF.MetaInfo(header.metainfo[1])
metainfo = VCF.MetaInfo(metainfo, tag="date", value="2022-08-17")
print(metainfo)
# add another line to the header like this:
push!(header, metainfo) 


writer = VCF.Writer(open("/home/annebusch/Documents/anne/PhD/test_vcf.vcf", "w"))
VCF.write(writer, header)
VCF.write(writer, record)
VCF.write(writer, records)

new_records = VCF.copy(records) 
VCF.write(writer, new_records)
VCF.write(writer, new_records)

close(writer)

# calling new_records.chrom returns the UnitRange{Int} in the new_records.data which describe the chromosome 
# this returns the name of the chromosome
print(VCF.chrom(new_records))
