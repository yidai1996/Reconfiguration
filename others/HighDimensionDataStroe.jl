using HDF5
# Writing data
A = zeros(48,15,4,4)
fid = h5open("filename.f5","w")
# create dataset
fid["Mydataset"] = A

# Reading data
fid2 = h5open("filename.f5","r")
B = read(fid2,"Mydataset")