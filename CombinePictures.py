import os
from PIL import Image
from fpdf import FPDF
 
# Import glob module to find all the files matching a pattern
import glob
 
PATH_WITH_FILES = "/home/jovyan/Outputs/Old/SmallerBatches"
image_extensions = ("*CalibrationCurve.png", "*CalibrationCurve.jpg", "*CalibrationCurve.gif")
output_file = "CalibrationCurves.pdf"

print("Starting print pictures to pdf program")

# This list will hold the images file names
images = []
 
# Build the image list by merging the glob results (a list of files)
# for each extension. We are taking images from current folder.
for extension in image_extensions:
    images.extend(glob.glob(PATH_WITH_FILES + "/" +  extension))
 

# Create instance of FPDF class
pdf=FPDF('P','in','letter')
# Add new page. Without this you cannot create the document.
pdf.add_page()
# Set font to Arial, 'B'old, 16 pts
pdf.set_font('Arial','B',16.0)
 
# Page header
pdf.cell(4.0,1.0,'Images in this folder')
pdf.ln(0.25)
 
# Smaller font for image captions
pdf.set_font('Arial','',10.0)
 
# Loop through the image list and position
# them with their caption one below the other
i = 1
for image in images:
 
    # Setting image width to half the page and
    # height to 1/4th of the page
    pdf.image(image,0,0, w=pdf.w/2.0, h=pdf.h/2.0)
    pdf.ln(0.15)
    
    # Image caption
    start = PATH_WITH_FILES + "/"
    end = '_'
    compound =  image[image.find(start)+len(start):image.rfind(end)]
    caption = "Figure " + i + ":  Compound: " +compound+" Moa: " 
    pdf.cell(3.0, 0.0, image)
    pdf.ln(0.25)
    i+=1
 
# output content into a file ('F')
pdf.output(output_file,'F')

print("Finished print pictures to pdf program")

# output_file = "CalibrationCurve.pdf"
# pdf = FPDF()
# PATH_WITH_FILES = "/home/jovyan/Outputs/Old/SmallerBatches"
# WELL_FILE_ENDING = "CalibrationCurve.png"
# w,h = 0,0
# i = 1

# for file_name in os.listdir(PATH_WITH_FILES):
#     if WELL_FILE_ENDING in file_name:
#     #    file_name = sdir + "IMG%.3d.png" % i
#         # if os.path.exists(file_name):
#             file_path =  PATH_WITH_FILES + "/" + file_name
#             if i == 1:
#                 cover = Image.open(file_path)
#                 w,h = cover.size
#                 pdf = FPDF(unit = "pt", format = [w,h])
#             image = file_path
#             pdf.add_page()
#             pdf.image(image,0,0,w,h)
#             i+=1
#         # else:
#         #     print("File not found:", file_name)
#             print("processed %d" % i)
# pdf.output(output_file, "F")
