-- This file will create 3D dipole interaction tensors
n=8
nx, ny, nz = n, n, n

lr = LongRange3D.new(nx,ny,nz)
ewald = DipoleEwald3D.new(nx,ny,nz)

-- fcc
sq = (1/2)^(1/2)
a = {0,sq,sq}
b = {sq,0,sq}
c = {sq,sq,0}
desc = "fcc"

ewald:setUnitCell(a, b, c)


ab = {"XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"}

mat = {}
for i=1,9 do
	print("Generating " .. ab[i] .. " elements")
	for z=0,nz-1 do
		for y=0,ny-1 do
			for x=0,nx-1 do
				-- Note the -1 factor here
 				local value = -1 * ewald:calculateTensorElement(ab[i], {x,y,z}) 
				lr:setMatrix(ab[i], {x,y,z}, value)
			end
		end
	end
end

filename = string.format("%dx%dx%d_%s.lua", nx, ny, nz, desc)
lr:saveTensors(filename)

print("Tensor saved to `" .. filename .. "' have a look at it, it's human readable")
