
# Nodos
def write_node_data(fname, nodes, data, nombre_datos=""):

    Ndata = len(nodes)
    fid = open(fname, "w")

    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$NodeData\n")
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0.\n") #eventualmente puede indicar el tiempo
    fid.write("3\n")
    fid.write("0\n")  #eventualmente puede indicar el paso tiempo
    fid.write("1\n")  # dimension de datos (1=escalares)
    fid.write(f"{Ndata}\n")  # numero de datos

    for i in range(Ndata):
        fid.write(f"{nodes[i]} {data[i]}\n")

    fid.write("$EndNodeData\n")

    fid.close()

    return

#Nodos
def write_node_data_2(fname, nodes, data1, data2, nombre_datos=""):

    Ndata = len(nodes)
    fid = open(fname, "w")

    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$NodeData\n")
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0.\n") #eventualmente puede indicar el tiempo
    fid.write("3\n")
    fid.write("0\n")  #eventualmente puede indicar el paso tiempo
    fid.write("3\n")  # dimension de datos (1=escalares)
    fid.write(f"{Ndata}\n")  # numero de datos

    for i in range(Ndata):
        fid.write(f"{nodes[i]} {data1[i]} {data2[i]} 0.0 \n")

    fid.write("$EndNodeData\n")

    fid.close()

    return

#elementos
def write_element_data(fname, elements, data, nombre_datos=""):

    Ndata = len(elements)
    fid = open(fname, "w")

    fid.write("$MeshFormat\n")
    fid.write("2.2 0 8\n")
    fid.write("$EndMeshFormat\n")
    fid.write("$ElementData\n")
    fid.write("1\n")
    fid.write(f"\"{nombre_datos}\"\n")
    fid.write("1\n")
    fid.write("0.\n")  # eventualmente puede indicar el tiempo
    fid.write("3\n")
    fid.write("0\n")  # eventualmente puede indicar el paso tiempo
    fid.write("1\n")  # dimension de datos (1=escalares)
    fid.write(f"{Ndata}\n")  # numero de datos

    for i in range(Ndata):
        fid.write(f"{elements[i]} {data[i]}\n")

    fid.write("$EndElementData\n")

    fid.close()


    return
