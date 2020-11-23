using Images
using Statistics
using ImageFiltering
using ImageTransformations
using Rotations
using CoordinateTransformations
using StaticArrays
using LinearAlgebra
using CoordinateTransformations
using BenchmarkTools
using OffsetArrays

include("Projective.jl")
function clp(c)
    return clamp(c,Gray(0.0),Gray(1.0))
end

xmin = 92
xmax = 2950
ymin = 90
ymax = 1810

i = 5

# Load image
img = load("source_images/$(lpad(i,3,'0')).jpg")
save("output_images/process/img.png",img)
h,w = size(img)

# Grayscale
img_gray = Gray.(img)
save("output_images/process/img_gray.png",img_gray)

# LoG filter
img_LoG = clp.(imfilter(img_gray, Kernel.LoG(3))*50)
save("output_images/process/img_LoG.png",img_LoG)

# Horizontal filter
n = 200
kernel_hrz = zeros(2n+1,2n+1)
kernel_hrz[n-1:n+1,:] .= 1/n
img_hrz = clp.(imfilter(img_LoG, kernel_hrz))
save("output_images/process/img_hrz.png",img_hrz)


# Vertical filter
n = 200
kernel_vrt = zeros(2n+1,2n+1)
kernel_vrt[:,n-1:n+1] .= 1/n
img_vrt = clp.(imfilter(img_LoG, kernel_vrt))
save("output_images/process/img_vrt.png",img_vrt)

# Multiply
img_mul = img_hrz.*img_vrt
save("output_images/process/img_mul.png",img_mul)

# Corner image
corner = 200
img_corner00 = img_mul[end-corner+1:end,1:corner]
img_corner10 = img_mul[end-corner+1:end,end-corner+1:end]
img_corner01 = img_mul[1:corner,1:corner]
img_corner11 = img_mul[1:corner,end-corner+1:end]

img_corner = mosaicview([img_corner01,img_corner00,img_corner11,img_corner10], npad=21, nrow=2, fillvalue=colorant"white")
save("output_images/process/img_corner.png",img_corner)

# Get coordinates of corners
corner_index00 = findmax(img_corner00)[2]+CartesianIndex(h-corner,0)
corner_index10 = findmax(img_corner10)[2]+CartesianIndex(h-corner,w-corner)
corner_index01 = findmax(img_corner01)[2]+CartesianIndex(0,0)
corner_index11 = findmax(img_corner11)[2]+CartesianIndex(0,w-corner)
corner_coordinates00 = [corner_index00[1],corner_index00[2]]
corner_coordinates10 = [corner_index10[1],corner_index10[2]]
corner_coordinates01 = [corner_index01[1],corner_index01[2]]
corner_coordinates11 = [corner_index11[1],corner_index11[2]]

# Generate marked image for debug
img_marked = RGB.(img_gray)
img_marked[[corner_index00
corner_index10
corner_index01
corner_index11]].=RGB(1,0,0)
img_marked[end-corner+1:end,1:corner]
img_marked[end-corner+1:end,end-corner-1:end]
img_marked[1:corner,1:corner]
img_marked[1:corner,end-corner+1:end]

# Transformation
mark00 = [ymax,xmin]
mark10 = [ymax,xmax]
mark01 = [ymin,xmin]
mark11 = [ymin,xmax]

p = Projective(mark00,mark10,mark01,mark11,corner_coordinates00,corner_coordinates10,corner_coordinates01,corner_coordinates11)
img_positioned = warp(img, p)
save("output_images/process/img_positioned.png",img_positioned)

# Crop
img_cropped = img_positioned[ymin-40:ymax+40,xmin-40:xmax+40]
save("output_images/process/img_cropped.png",img_cropped)

# Resize
img_resized = imresize(img_cropped, ratio=1/5)
save("output_images/process/img_resized.png",img_resized)

# Color alignment
img_target = load("source_images/001.jpg")
meancolor = mean(img_target)
img_coloraligned = img_resized.-mean(img_resized).+meancolor
save("output_images/process/img_coloraligned.png",img_coloraligned)
