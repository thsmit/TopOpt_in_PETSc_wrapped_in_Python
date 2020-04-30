#!/usr/bin/python3

# Modified by: Thijs Smit, May 2020

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

import struct as st
import binascii as ba

FOUT = "vtutest.vtu"

def main():
	try:
		fout = open(FOUT,'wb')
	except:
		exit("Couldn't create file...")

	writeHeader(fout,4,1)
	lP = [0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 1, 1]
	writePoints(fout,lP)
	lConn =[0, 1, 2, 3]
	lOffset = [4]
	lTypes = [10]
	writeCells(fout,lConn,lOffset,lTypes)
	lldata = []
	lldata.append(["Scalars",1,2,-1,0.5])
	lldata.append(["OtherSet",1.0,2,3,127.3])
	writeScalarPointData(fout,lldata)
	lldata = []
	lldata.append(["Scalars",1.0])
	lldata.append(["Scalars2",27.0])
	writeScalarCellData(fout,lldata)
	writeFooter(fout)
	fout.close()

#nP = number of points, nC = number of cells
def writeHeader(fout,nP,nC):
	fout.write("<?xml version=\"1.0\"?>\n".encode('ascii'))
	fout.write("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" header_type=\"UInt64\" byte_order=\"LittleEndian\">\n".encode('ascii'))
	fout.write("<UnstructuredGrid>\n".encode('ascii'))
	fout.write("\t<Piece NumberOfPoints=\"".encode('ascii') + st.pack("i", nP) + "\" NumberOfCells=\"".encode('ascii') + st.pack("i", nC) + "\">\n".encode('ascii'))


def writeFooter(fout):
	fout.write("\t</Piece>\n".encode('ascii'))
	fout.write("</UnstructuredGrid>\n".encode('ascii'))
	fout.write("</VTKFile>".encode('ascii'))


def writeRawPoints(fout,data):
	fout.write("\t\t<Points>\n".encode('ascii'))
	fout.write("\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n".encode('ascii'))
	fout.write(writeBin64(data))
	fout.write("\t\t\t</DataArray>\n".encode('ascii'))
	fout.write("\t\t</Points>\n".encode('ascii'))


#lP = list of points
#3D coordinates assumed...
def writePoints(fout,lP):
	fout.write("\t\t<Points>\n".encode('ascii'))
	fout.write("\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n".encode('ascii'))
	data = ""
	for number in lP:
		#Append number to binary datalist
		data += st.pack('f',number)
	fout.write(writeBin64(data))
	fout.write("\t\t\t</DataArray>\n".encode('ascii'))
	fout.write("\t\t</Points>\n".encode('ascii'))

#lConn = connectivity list
#lOffset = Offset list
def writeCells(fout,lConn,lOffset,lTypes):
	fout.write("\t\t<Cells>\n".encode('ascii'))
	fout.write("\t\t\t<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"binary\">\n".encode('ascii'))
	data = ""
	for number in lConn:
		data += st.pack('I',number)
	fout.write(writeBin64(data))
	fout.write("\t\t\t</DataArray>\n".encode('ascii'))

	fout.write("\t\t\t<DataArray type=\"UInt32\" Name=\"offsets\" format=\"binary\">\n".encode('ascii'))
	data = ""
	for number in lOffset:
		data += st.pack('I',number)
	fout.write(writeBin64(data))
	fout.write("\t\t\t</DataArray>\n".encode('ascii'))

	#Data type could be UInt8 instead, but UInt32 keeps it simple for now...
	fout.write("\t\t\t<DataArray type=\"UInt32\" Name=\"types\" format=\"binary\">\n".encode('ascii'))
	data = ""
	for number in lTypes:
		data += st.pack('I',number)
	fout.write(writeBin64(data))
	fout.write("\t\t\t</DataArray>\n".encode('ascii'))
	fout.write("\t\t</Cells>".encode('ascii'))



def writeRawCellsConn(fout,data):
	fout.write("\t\t<Cells>\n".encode('ascii'))
	fout.write("\t\t\t<DataArray type=\"UInt64\" Name=\"connectivity\" format=\"binary\">\n".encode('ascii'))
	fout.write(writeBin64(data))
	fout.write("\t\t\t</DataArray>\n".encode('ascii'))

def writeRawCellsOffset(fout,data):
	fout.write("\t\t\t<DataArray type=\"UInt64\" Name=\"offsets\" format=\"binary\">\n".encode('ascii'))
	fout.write(writeBin64(data))
	fout.write("\t\t\t</DataArray>\n".encode('ascii'))

def writeRawCellsType(fout,data):
	fout.write("\t\t\t<DataArray type=\"UInt64\" Name=\"types\" format=\"binary\">\n".encode('ascii'))
	fout.write(writeBin64(data))
	fout.write("\t\t\t</DataArray>\n".encode('ascii'))
	fout.write("\t\t</Cells>".encode('ascii'))

# Returns Int containing data size + data
# everything in base64 encoding
def writeBin64(bindata):
	datalen = st.pack('Q',len(bindata))
	# Note that \n is removed from datalen
	# - which btw is a ridiculous design choice
	# as it makes the data size grow by 33% (4/3-1)
	return ba.b2a_base64(datalen)[0:-1]+ba.b2a_base64(bindata.encode('ascii'))


##lldata = list of list of data (where 1st entry is set name)
def writeScalarPointData(fout,lldata):
	fout.write("\t\t<PointData Scalars=\"scalars\">\n".encode('ascii'))

	for ldata in lldata:
		setname = ldata[0]
		del ldata [0]
		data = ""
		for number in ldata:
			data += st.pack('f',number)
		fout.write("\t\t\t<DataArray type=\"Float32\" Name=\"".encode('ascii') + setname + "\" format=\"binary\">\n".encode('ascii'))
		fout.write(writeBin64(data))
		fout.write("\t\t\t</DataArray>\n".encode('ascii'))
	
	fout.write("\t\t</PointData>\n".encode('ascii'))


def writeRawScalarPointData(fout,ldata,setnames):
	fout.write("\t\t<PointData Scalars=\"scalars\">\n".encode('ascii'))
	j=-1
	for data in ldata:
		j=j+1
		fout.write("\t\t\t<DataArray type=\"Float32\" Name=\"".encode('ascii') + setnames[j].encode('ascii') + "\" format=\"binary\">\n".encode('ascii'))
		fout.write(writeBin64(data))
		fout.write("\t\t\t</DataArray>\n".encode('ascii'))
	
	fout.write("\t\t</PointData>\n".encode('ascii'))


def writeScalarCellData(fout,lldata):
	fout.write("\t\t<CellData Scalars=\"scalars\">\n".encode('ascii'))

	for ldata in lldata:
		setname = ldata[0]
		del ldata [0]
		data = ""
		for number in ldata:
			data += st.pack('f',number)
		fout.write("\t\t\t<DataArray type=\"Float32\" Name=\"".encode('ascii') + setname.encode('ascii') + "\" format=\"binary\">\n".encode('ascii'))
		fout.write(writeBin64(data))
		fout.write("\t\t\t</DataArray>\n".encode('ascii'))

	fout.write("\t\t</CellData>\n".encode('ascii'))


def writeRawScalarCellData(fout,ldata,setnames):
	fout.write("\t\t<CellData Scalars=\"scalars\">\n".encode('ascii'))

	j=-1
	for data in ldata:
		j=j+1
		fout.write("\t\t\t<DataArray type=\"Float32\" Name=\"".encode('ascii') + setnames[j].encode('ascii') + "\" format=\"binary\">\n".encode('ascii'))
		fout.write(writeBin64(data))
		fout.write("\t\t\t</DataArray>\n".encode('ascii'))

	fout.write("\t\t</CellData>\n".encode('ascii'))


if __name__ == "__main__":
	main()
