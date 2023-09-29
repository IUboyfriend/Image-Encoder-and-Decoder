import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mping
import os
import math
def readfile():
    while(1):#The program will loop until a valid input is received. 
        filename=input("Please input the file name of an image:\n")
        try:#Use try...except... to handle the exceptions in case that the file does not exist.
            img=mping.imread(filename)
        except IOError:
            print("<{}> does not exist!".format(filename))
        else:
            return filename
def readd():
    while(1):#The program will loop until a valid input is received. 
        d=input("Please input the value of d:\n")
        if(d==''):#If null string is input, return the default value 4 of d.
            return 4
        if(int(d)>1):#check whether d is an integer power of 2
            i=1
            while(i<=int(d)):
                if(i==int(d)):
                    return int(d)
                i*=2
        print("d must be an integer power of 2!\n")
def codeword(a,g0,g1):
    minn=(a[0,0]-g0)**2+(a[0,1]-g0)**2+(a[1,0]-g1)**2+(a[1,1]-g1)**2#Calculate the distance between the subblock and codewords C0, then give this value to the variable minn to record the current min value.
    code=0#Initial the codeword to be chosen to codeword C0
    count1=(a[0,0]-g1)**2+(a[0,1]-g1)**2+(a[1,0]-g0)**2+(a[1,1]-g0)**2#Calculate the distance between the subblock and codewords C1.
    if(count1<minn):#If this distance is less than minn, change the minn and code accordingly.
        minn=count1
        code=1
    count2=(a[0,0]-g0)**2+(a[0,1]-g1)**2+(a[1,0]-g0)**2+(a[1,1]-g1)**2#Calculate the distance between the subblock and codewords C2.
    if(count2<minn):#If this distance is less than minn, change the minn and code accordingly.
        minn=count2
        code=2
    count3=(a[0,0]-g1)**2+(a[0,1]-g0)**2+(a[1,0]-g1)**2+(a[1,1]-g0)**2#Calculate the distance between the subblock and codewords C3.
    if(count3<minn):#If this distance is less than minn, change the minn and code accordingly.
        minn=count3
        code=3
    return code
        
def BVQCencode(in_image_filename,out_encoding_result_filename,d):
    img=mping.imread(in_image_filename)#read the image
    plt.imshow(img,cmap='gray')
    plt.show()#show the image
    X=np.array(img)*255#Normalise the pixel values to (0,255).
    rows_of_blocks=int(X.shape[0]/d)#Calculate the number of rows for blocks
    cols_of_blocks=int(X.shape[1]/d)#Calculate the number of columns for blocks
    meanplane=np.zeros([rows_of_blocks,cols_of_blocks],dtype="uint8")#create an array to store the mean of each block
    stdplane=np.zeros([rows_of_blocks,cols_of_blocks],dtype="uint8")#create an array to store the standard deviation of each block
    codeplane=np.zeros([int(X.shape[0]/2),int(X.shape[1]/2)],"uint8")#create an array to store the code words 
    for i in range(0,X.shape[0],d):
        for j in range(0,X.shape[1],d):#enumerate each block
            block=X[i:i+d,j:j+d]#get the block from the normalised image
            miu=np.mean(block)#calculate the mean of the current block
            std=np.std(block)#calculate the standard deviation of the current block
            meanplane[int(i/d),int(j/d)]=np.uint8(round(miu))#store the mean value of this block to the appropriate position of the meanplane
            stdplane[int(i/d),int(j/d)]=np.uint8(round(std))#store the standard deviation value of this block to the appropriate position of the stdplane
            g0=max(0,miu-std)#compute g0
            g1=min(255,miu+std)#compute g1
            for m in range(0,d,2):
                for n in range(0,d,2): #enumerate each subblock in the current block
                    code=codeword(X[i+m:i+m+2,j+n:j+n+2],g0,g1)#get the code for the current subblock
                    codeplane[int((i+m)/2),int((j+n)/2)]=np.uint8(code)##store the code of the subblock to the appropriate position of the codeplane
    file=open(out_encoding_result_filename,"wb")
    header=np.zeros([6],dtype='uint8')#write the header
    header[0]=np.uint8(6)
    header[1]=np.uint8(d)
    header[2]=np.uint8(cols_of_blocks%256)
    header[3]=np.uint8(cols_of_blocks/256)
    header[4]=np.uint8(rows_of_blocks%256)
    header[5]=np.uint8(rows_of_blocks/256)
    for byte in header:
        file.write(byte)
    for i in range(0,rows_of_blocks):
        for j in range(0,cols_of_blocks):#enumerate each block
            file.write(meanplane[i,j])#write the mean of this block into the file
            file.write(stdplane[i,j])#write the standard deviation of this block into the file
            cnt=0#used to count from 0 to 4 so that we can know when a byte will be generated
            now=0#used to store the current value of the byte
            for m in range(0,d,2):
                for n in range(0,d,2):#enumerate each subblock
                    cnt+=1
                    now=now*4+codeplane[int((i*d+m)/2),int((j*d+n)/2)]#change the current value of the byte
                    if cnt==4:#if 4 operations have been done, write the byte into the file and initialise the cnt and now
                        file.write(np.uint8(now))
                        cnt=0
                        now=0
            if cnt!=0:#if finally, the cnt is not 0, it means that we have to add K bits to form a complete byte,
                now*=(4**(4-cnt))
                file.write(np.uint8(now))
    file.close()
    return ({'M':meanplane,'Sd':stdplane,'Idx':codeplane})

def BVQCdecode(in_encoding_result_filename, out_reconstructed_image_filename):
    file=open(in_encoding_result_filename,"rb")
    header_len=file.read(1)[0]#read the header
    d=file.read(1)[0]
    no_of_block_cols=file.read(1)[0]+file.read(1)[0]*256
    no_of_block_rows=file.read(1)[0]+file.read(1)[0]*256
    OImg=np.zeros([no_of_block_rows*d,no_of_block_cols*d])#create an array to store the pixel values after decoding
    n=np.ceil(d/2*d/2*2/8)+2#n represents the number of bytes we have to read for each block 
    for i in range(0,no_of_block_rows):
        for j in range(0,no_of_block_cols):#enumerate each block
            bfr=file.read(np.uint(n))#read one content for a block
            mean=bfr[0]#read the mean
            std=bfr[1]#read the standard deviation
            g0=max(0,mean-std)#calculate g0
            g1=min(255,mean+std)#calculate g1
            cnt=0#cnt is used to store how many indexes have been read
            for byte in bfr[2:]:
                temp=np.uint8(byte)
                for m in range(0,4):#deal with the four indexes in a byte
                    current=int (temp/(64/(4**m)))#get the (m+1)th index in a byte
                    temp-=64/(4**m)*current#modify twmp accordingly
                    now_row=int(cnt/(d/2))#based on the number of indexes that have been read, we can define the index is in which row and column in the codeplane
                    now_col=int(cnt%(d/2))
                    #Based on the index, rebulid the block with according values in the code book
                    if(current==0):
                        OImg[i*d+2*now_row,j*d+2*now_col]=g0
                        OImg[i*d+2*now_row,j*d+2*now_col+1]=g0
                        OImg[i*d+2*now_row+1,j*d+2*now_col]=g1
                        OImg[i*d+2*now_row+1,j*d+2*now_col+1]=g1
                    if(current==1):
                        OImg[i*d+2*now_row,j*d+2*now_col]=g1
                        OImg[i*d+2*now_row,j*d+2*now_col+1]=g1
                        OImg[i*d+2*now_row+1,j*d+2*now_col]=g0
                        OImg[i*d+2*now_row+1,j*d+2*now_col+1]=g0
                    if(current==2):
                        OImg[i*d+2*now_row,j*d+2*now_col]=g0
                        OImg[i*d+2*now_row,j*d+2*now_col+1]=g1
                        OImg[i*d+2*now_row+1,j*d+2*now_col]=g0
                        OImg[i*d+2*now_row+1,j*d+2*now_col+1]=g1
                    if(current==3):
                        OImg[i*d+2*now_row,j*d+2*now_col]=g1
                        OImg[i*d+2*now_row,j*d+2*now_col+1]=g0
                        OImg[i*d+2*now_row+1,j*d+2*now_col]=g1
                        OImg[i*d+2*now_row+1,j*d+2*now_col+1]=g0
                    cnt+=1
                    if cnt==int(d/2*d/2):#This part is used to handle the case that stuffing bits exists.
                        break
    file.close()
    plt.imshow(OImg,cmap='gray')
    plt.show()#show the reconstructed image
    plt.imsave(out_reconstructed_image_filename,OImg,cmap='gray')  #save the image 
    return OImg
def evaluate(original,product):
    size=X.shape[0]*X.shape[1]
    ans=0
    for i in range(0,X.shape[0]):
        for j in range(0,X.shape[1]):
            ans+=(original[i][j]-product[i][j])**2#for each pixel in original image and the reconstructed image,calculate the error.
    mse=ans/size
    PPSNR=10*math.log(255*255/mse,10)
    return mse,PPSNR
    
    
    
in_file=readfile()#Read the file.
out_file="image_encoded.out"#The encoded image is stored in this file.
d=readd()#Store d.
result=BVQCencode(in_file,out_file,d)#  Do the encode and the returned dictionary is stored in variable result. 
print(result) 
file_in="image_encoded.out"#Change the input file and output file name to do the decode.
file_out="image_decoded.png"
OImg=BVQCdecode(file_in,file_out) # OImg is an array store pixel values of the decoded image. 
img=mping.imread(in_file)#This in_file stored the original pixel values of the image which is bounded in (0,1).
X=np.array(img)*255#Normalise the pixel values to (0,255) to calculate the error.
MSE,PPSNR=evaluate(X,OImg)    #Calculate the error.
print("MSE={}".format(MSE)) 
print("PPSNR={}".format(PPSNR))
