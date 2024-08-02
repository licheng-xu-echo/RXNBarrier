from ase.io import read as ase_read
def read_xyz(file_path):
    with open(file_path,'r') as fr:
        lines = fr.readlines()
    atom_num = int(lines[0].strip())
    symbol_coords_lines = lines[2:2+atom_num]
    symbols = [line.split()[0] for line in symbol_coords_lines]
    coords = [list(map(float,line.split()[1:4])) for line in symbol_coords_lines]
    return symbols,coords