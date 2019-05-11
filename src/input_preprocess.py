import math
filename="D:\\test\\Bridge-4.inp"
f=open(filename,mode="r")
fl=f.readlines()
row_num_f=len(fl)
i=0
keycase=1
part_list=[]
part_id=0
node_list=[]
element_list=[]
instance_list=[]
tie_nset_list=[]
Boundary_nset_list=[]
# the node data of the part
part_node_dict={"test":0}
# the element type of the part
part_eletype_dict={"test":0}
# the element data of the part
part_element_dict={"test":0}
# the node set name of the part
part_nset_dict={"test":0}
# the element set name of the part
part_elset_dict={"test":0}
# the section type of the part
part_section_name_dict={"test":0}
# the node id of the node set in the part
part_nset_dict={"test":0}
# the element id of the element set in the part
part_elset_dict={"test":0}
# the material name of the element in the part
part_section_material_dict={"test":0}
# the parameter of the material of the element in the part
part_material_para_dict={"test":0}
# the part name of the instance
instance_part_dict={"test":0}
# the node data of the instance
instance_node_dict={"test":0}
# the element data of the instance
instance_element_dict={"test":0}
# the target instance of the nset
nset_instance_dict={"test":0}
# the target node id of the nset in the instance 
nset_node_dict={"test":0}
# the nset to form the surface
surface_nset_dict={"test":0}
# the density of the material
material_density_dict={"test":0}
# the E of the material
material_E_dict={"test":0}
# the poisson value of the material
material_poisson_dict={"test":0}
# the start degree of the boundary node defined
nset_start_degree1_dict={"test":0}
nset_start_degree2_dict={"test":0}
nset_start_degree3_dict={"test":0}
# the end degree of the boundary node defined
nset_end_degree1_dict={"test":0}
nset_end_degree2_dict={"test":0}
nset_end_degree3_dict={"test":0}
# the replace dictionary for the tied points
node_replace_dict={"test":0}
def myrange(n):
	ss=0
	list_range=[]
	while ss<n:
		list_range.append(ss+1)
		ss=ss+1
	return list_range
def get_list_replace(n):
	list_n=myrange(n)
	for kk in set(node_replace_dict.keys()):
		if kk=="test":
			pass
		else:
			ll=int(kk)
			while ll<n:
				list_n[ll]=list_n[ll]-1
				ll=ll+1
	print(list_n)
	return list_n
while i<row_num_f:
	row=fl[i]
	row=row.rstrip("\n")
	row=row.replace(" ","")
	if row=="":
		keyold=keycase
		keycase=0
	else:
		# get the control word
		len_row=len(row)
		if len_row<1:
			control_ch1=""
			control_ch2=""
		elif len_row==1:
			control_ch1=row[0]
			control_ch2=""
		else:
			control_ch1=row[0]
			control_ch2=row[1]
	# keycase 0: Notes line
	if control_ch1=="*":
		if control_ch2=="*":
			keyold=keycase
			keycase=0
		else:
			pass
	else:
		pass
	# notes line
	if keycase==0:
		keycase=keyold
	# control line
	elif keycase==1:
		if control_ch1=="*":
			key_str=row[1:]
			key_list=key_str.split(",")
			keyword=key_list[0]
			if keyword=="Heading":
				pass
			elif keyword=="Preprint":
				pass
			elif keyword=="Part":
				keycase=2
				part_name=key_list[1][5:]
				part_list.append(part_name)
				part_id=part_id+1
				node_list.append([])
				element_list.append([])
				part_element_dict[part_name]=[]
				part_node_dict[part_name]=[]
			elif keyword=="Assembly":
				keycase=12
			elif keyword=="Material":
				control_str1=key_list[1].split("=")
				if control_str1[0]=="name":
					material_name=control_str1[1]
				else:
					print("error:new case in the name of material")
				keycase=18
			elif keyword=="Step":
				keycase=21
	# read the part data
	elif keycase==2:
		if control_ch1=="*":
			key_str=row[1:]
			key_list=key_str.split(",")
			keyword=key_list[0]
			if keyword=="EndPart":
				keycase=1
			elif keyword=="Node":
				keycase=3
			elif keyword=="Element":
				type_name=key_list[1].split("=")[1]
				part_eletype_dict[part_name]=type_name
				keycase=4
			elif keyword=="Nset":
				control_str=key_list[1].split("=")
				if control_str[0]=="nset":
					part_nset_dict[part_name]=control_str[1]
					nset_name=control_str[1]
				else:
					print("error:occur other instance")
					print(i)
				generate_buel=key_list[len(key_list)-1]
				if generate_buel=="generate":
					keycase=5
				else:
					keycase=6
			elif keyword=="Elset":
				control_str=key_list[1].split("=")
				if control_str[0]=="elset":
					part_elset_dict[part_name]=control_str[1]
					elset_name=control_str[1]
				else:
					print("error:occur other instance")
					print(i)
				generate_buel=key_list[len(key_list)-1]
				if generate_buel=="generate":
					keycase=7
				else:
					keycase=8
			elif keyword=="SolidSection":
				part_section_name_dict[part_name]="SolidSection"
				control_str2=key_list[2].split("=")
				if control_str2[0]=="material":
					part_section_material_dict[part_name]=control_str2[1]
				else:
					print("error: new case")
					print(i)
				keycase=9
			elif keyword=="ShellSection":
				part_section_name_dict[part_name]="ShellSection"
				control_str2=key_list[2].split("=")
				if control_str2[0]=="material":
					part_section_material_dict[part_name]=control_str2[1]
				else:
					print("error: new case")
					print(i)
				keycase=10
			elif keyword=="BeamSection":
				control_str3=key_list[len(key_list)-1].split("=")
				if control_str3[0]=="section":
					beam_section_type=control_str3[1]
				part_section_name_dict[part_name]="BeamSection"+"_"+beam_section_type
				control_str2=key_list[2].split("=")
				if control_str2[0]=="material":
					part_section_material_dict[part_name]=control_str2[1]
				else:
					print("error: new case")
					print(i)
				if beam_section_type=="BOX":
					keycase=11
					part_material_para_dict[part_name]=[]
		else:
			print("Error: false input after the part")
			print(i)
	# read the node data
	elif keycase==3:
		if control_ch1=="*":
			i=i-1
			keycase=2
		else:
			node_list[part_id-1].append(row.split(","))
			part_node_dict[part_name].append(row.split(","))
	elif keycase==4:
		if control_ch1=="*":
			i=i-1
			keycase=2
		else:
			element_list[part_id-1].append(row.split(","))
			part_element_dict[part_name].append(row.split(","))
	elif keycase==5:
		if control_ch1=="*":
			i=i-1
			keycase=2
		else:
			node_list_temp=row.split(",")
			ik=int(node_list_temp[0])
			jk=int(node_list_temp[1])
			pace=int(node_list_temp[2])
			nodeid_list=[]
			k=ik
			while k<jk+1:
				nodeid_list.append(k)
				k=k+pace
			part_nset_dict[part_name]=nodeid_list
	elif keycase==6:
		if control_ch1=="*":
			i=i-1
			keycase=2
		else:
			node_list_temp=row.split(",")
			nodeid_list=[]
			for id in node_list_temp:
				nodeid_list.append(int(id))
			part_nset_dict[part_name]=nodeid_list
	elif keycase==7:
		if control_ch1=="*":
			i=i-1
			keycase=2
		else:
			element_list_temp=row.split(",")
			ik=int(element_list_temp[0])
			jk=int(element_list_temp[1])
			pace=int(element_list_temp[2])
			elementid_list=[]
			k=ik
			while k<jk+1:
				elementid_list.append(k)
				k=k+pace
			part_elset_dict[part_name]=elementid_list
	elif keycase==8:
		if control_ch1=="*":
			i=i-1
			keycase=2
		else:
			element_list_temp=row.split(",")
			elementid_list=[]
			if element_list_temp[len(element_list_temp)-1]=="":
				element_list_temp=element_list_temp[:len(element_list_temp)-1]
			for id in element_list_temp:
				elementid_list.append(int(id))
			part_elset_dict[part_name]=elementid_list
	elif keycase==9:
		if control_ch1=="*":
			i=i-1
			keycase=2
		else:
			material_list=row.split(",")
			if material_list[len(material_list)-1]=="":
				material_list=material_list[:len(material_list)-1]
			part_material_para_dict[part_name]=material_list
	elif keycase==10:
		if control_ch1=="*":
			i=i-1
			keycase=2
		else:
			material_list=row.split(",")
			if material_list[len(material_list)-1]=="":
				material_list=material_list[:len(material_list)-1]
			part_material_para_dict[part_name]=material_list
	elif keycase==11:
		if control_ch1=="*":
			i=i-1
			keycase=2
		else:
			material_list=row.split(",")
			if material_list[len(material_list)-1]=="":
				material_list=material_list[:len(material_list)-1]
			part_material_para_dict[part_name].append(material_list)
	elif keycase==12:
		if control_ch1=="*":
			key_str=row[1:]
			key_list=key_str.split(",")
			keyword=key_list[0]
			if keyword=="Instance":
				control_str1=key_list[1].split("=")
				if control_str1[0]=="name":
					instance_name=control_str1[1]
				else:
					print("error: new case in instance name")
				instance_list.append(instance_name)
				control_str2=key_list[2].split("=")
				if control_str2[0]=="part":
					part_name=control_str2[1]
				else:
					print("error: new case in instance name")
				instance_part_dict[instance_name]=part_name
				instance_transform_list=[]
				instance_node_dict[instance_name]=[]
				keycase=13
			elif keyword=="Nset":
				keycase=14
				control_str1=key_list[1].split("=")
				if control_str1[0]=="nset":
					nset_name=control_str1[1]
				else:
					print("error: new case in the nset name of the Assembly")
				control_str2=key_list[2].split("=")
				if control_str2[0]=="instance":
					instance_name=control_str2[1]
				else:
					print("error: new case int the target instance name of the nset(Assembly)")
				if nset_name in nset_instance_dict.keys():
					if isinstance(nset_instance_dict[nset_name],list):
						nset_instance_dict[nset_name].append(instance_name)
					else:
						nset_instance_dict[nset_name]=[nset_instance_dict[nset_name]]
				else:
					nset_instance_dict[nset_name]=[instance_name]
				help_nset=[]
			elif keyword=="Elset":
				keycase=15
			elif keyword=="Surface":
				control_str1=key_list[1].split("=")
				if control_str1[0]=="type":
					if control_str1[1]=="NODE":
						pass
					else:
						print("error:the Surface is not defined from node")
						print(i)
				else:
					print("new case in Surface type")
				control_str2=key_list[2].split("=")
				if control_str2[0]=="name":
					surface_name=control_str2[1]
				else:
					print("error new case in surface name")
				control_str3=key_list[3]
				if control_str3=="internal":
					pass
				else:
					print("error: the internal buel not on")
				keycase=16
			elif keyword=="Tie":
				keycase=17
			elif keyword=="EndAssembly":
				keycase=1
		else:
			print("error:false input after Assembly")
			print(i)
	elif keycase==13:
		if control_ch1=="*":
			key_str=row[1:]
			key_list=key_str.split(",")
			keyword=key_list[0]
			if keyword=="EndInstance":
				keycase=12
			else:
				print("error:new case in end instance")
			part_node_list_all=part_node_dict[part_name]
			if len(instance_transform_list)==0:
				for part_node_list in part_node_list_all:
					instance_node_id=int(float(part_node_list[0]))
					instance_node_x=float(part_node_list[1])
					instance_node_y=float(part_node_list[2])
					instance_node_z=float(part_node_list[3])
					instance_node_list=[instance_node_id,instance_node_x,instance_node_y,instance_node_z,[],0]
					instance_node_dict[instance_name].append(instance_node_list)
			elif len(instance_transform_list)==1:
				translate_list=instance_transform_list[0]
				translate_x=float(translate_list[0])
				translate_y=float(translate_list[1])
				translate_z=float(translate_list[2])
				for part_node_list in part_node_list_all:
					part_node_x=float(part_node_list[1])
					part_node_y=float(part_node_list[2])
					part_node_z=float(part_node_list[3])
					instance_node_id=int(float(part_node_list[0]))
					instance_node_x=part_node_x+translate_x
					instance_node_y=part_node_y+translate_y
					instance_node_z=part_node_z+translate_z
					instance_node_list=[instance_node_id,instance_node_x,instance_node_y,instance_node_z,[],0]
					instance_node_dict[instance_name].append(instance_node_list)
			elif len(instance_transform_list)==2:
				translate_list=instance_transform_list[0]
				translate_x=float(translate_list[0])
				translate_y=float(translate_list[1])
				translate_z=float(translate_list[2])
				rotation_list=instance_transform_list[1]
				if len(rotation_list)!=7:
						print("error:new case in rotation_list")
						print(rotation_list)
				a_x=float(rotation_list[0])
				a_y=float(rotation_list[1])
				a_z=float(rotation_list[2])
				b_x=float(rotation_list[3])
				b_y=float(rotation_list[4])
				b_z=float(rotation_list[5])
				theta=float(rotation_list[6])*math.pi/180.0
				for part_node_list in part_node_list_all:
					part_node_x=float(part_node_list[1])
					part_node_y=float(part_node_list[2])
					part_node_z=float(part_node_list[3])
					if len(translate_list)!=3:
						print("error:new case in translate_list")
					else:
						pass
					instance_node_x_translate=part_node_x+translate_x
					instance_node_y_translate=part_node_y+translate_y
					instance_node_z_translate=part_node_z+translate_z
					instance_node_id=int(float(part_node_list[0]))
					instance_node_x=a_x + ((math.sin(theta)*(a_z - b_z))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2)**(1/2) - ((math.cos(theta) - 1)*(a_x - b_x)*(a_y - b_y))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2))*(instance_node_y_translate - a_y) - ((math.sin(theta)*(a_y - b_y))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2)**(1/2) + ((math.cos(theta) - 1)*(a_x - b_x)*(a_z - b_z))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2))*(instance_node_z_translate - a_z) + (instance_node_x_translate - a_x)*(math.cos(theta) - ((math.cos(theta) - 1)*(a_x - b_x)**2)/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2))
					instance_node_y=a_y - ((math.sin(theta)*(a_z - b_z))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2)**(1/2) + ((math.cos(theta) - 1)*(a_x - b_x)*(a_y - b_y))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2))*(instance_node_x_translate - a_x) + ((math.sin(theta)*(a_x - b_x))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2)**(1/2) - ((math.cos(theta) - 1)*(a_y - b_y)*(a_z - b_z))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2))*(instance_node_z_translate - a_z) + (instance_node_y_translate - a_y)*(math.cos(theta) - ((math.cos(theta) - 1)*(a_y - b_y)**2)/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2))
					instance_node_z=a_z + ((math.sin(theta)*(a_y - b_y))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2)**(1/2) - ((math.cos(theta) - 1)*(a_x - b_x)*(a_z - b_z))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2))*(instance_node_x_translate - a_x) - ((math.sin(theta)*(a_x - b_x))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2)**(1/2) + ((math.cos(theta) - 1)*(a_y - b_y)*(a_z - b_z))/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2))*(instance_node_y_translate - a_y) + (instance_node_z_translate - a_z)*(math.cos(theta) - ((math.cos(theta) - 1)*(a_z - b_z)**2)/((a_x - b_x)**2 + (a_y - b_y)**2 + (a_z - b_z)**2))
					instance_node_list=[instance_node_id,instance_node_x,instance_node_y,instance_node_z,[],0]
					instance_node_dict[instance_name].append(instance_node_list)
			else:
				print("error:new case in instance data")
				print(i)
				print(instance_transform_list)
		else:
			transform_data=row.split(",")
			transform_data_list=[]
			for data_transform in transform_data:
				transform_data_list.append(float(data_transform))
			instance_transform_list.append(transform_data_list)
	elif keycase==14:
		if control_ch1=="*":
			i=i-1
			keycase=12
			if nset_name in nset_node_dict.keys():
					if isinstance(nset_node_dict[nset_name][0],list):
						nset_node_dict[nset_name].append(help_nset)
					else:
						print("error in keycase 14")
			else:
				nset_node_dict[nset_name]=[help_nset]
		else:
			nset_node_list_data=row.split(",")
			if nset_node_list_data[len(nset_node_list_data)-1]=="":
				nset_node_list_data=nset_node_list_data[:(len(nset_node_list_data)-1)]
			else:
				pass
			nset_node_id_list=[]
			for nset_node_list_data_element in nset_node_list_data:
				nset_node_id_list.append(int(nset_node_list_data_element))
			help_nset=help_nset+nset_node_id_list
	elif keycase==15:
		if control_ch1=="*":
			i=i-1
			keycase=12
		else:
			pass
	elif keycase==16:
		if control_ch1=="*":
			i=i-1
			keycase=12
		else:
			if float(row.split(",")[1])==1.0:
				surface_nset_dict[surface_name]=row.split(",")[0]
			else:
				print("the default area of the surface is not unit")
	elif keycase==17:
		if control_ch1=="*":
			i=i-1
			keycase=12
		else:
			tie_nset_list.append(row.split(","))
	elif keycase==18:
		if control_ch1=="*":
			key_str=row[1:]
			key_list=key_str.split(",")
			keyword=key_list[0]
			if keyword=="Density":
				keycase=19
			elif keyword=="Elastic":
				keycase=20
			elif keyword=="Material":
				i=i-1
				keycase=1
			elif keyword=="Step":
				i=i-1
				keycase=1
			else:
				print("error:wrong input after the material")
		else:
			print("error:wrong input after the material")
			print(row)
			print(i)
	elif keycase==19:
		if control_ch1=="*":
			i=i-1
			keycase=18
		else:
			material_density_dict[material_name]=float(row.split(",")[0])
	elif keycase==20:
		if control_ch1=="*":
			i=i-1
			keycase=18
		else:
			material_E_dict[material_name]=float(row.split(",")[0])	
			material_poisson_dict[material_name]=float(row.split(",")[1])	
	elif keycase==21:
		if control_ch1=="*":
			key_str=row[1:]
			key_list=key_str.split(",")
			keyword=key_list[0]
			if keyword=="Static":
				keycase=22
			elif keyword=="Boundary":
				keycase=23
			elif keyword=="Dload":
				keycase=24
			elif keyword=="Restart":
				keycase=24
			elif keyword=="Output":
				keycase=24
			elif keyword=="EndStep":
				pass
	elif keycase==22:
		if control_ch1=="*":
			i=i-1
			keycase=21
		else:
			pass
	elif keycase==23:
		if control_ch1=="*":
			i=i-1
			keycase=21
		else:
			Boundary_list=row.split(",")
			Boundary_nset=Boundary_list[0]
			Boundary_nset_list.append(Boundary_nset)
			if Boundary_list[1]=="1":
				nset_start_degree1_dict[Boundary_nset]=Boundary_list[1]
				nset_end_degree1_dict[Boundary_nset]=Boundary_list[2]
			elif Boundary_list[1]=="2":
				nset_start_degree2_dict[Boundary_nset]=Boundary_list[1]
				nset_end_degree2_dict[Boundary_nset]=Boundary_list[2]
			elif Boundary_list[1]=="3":
				nset_start_degree3_dict[Boundary_nset]=Boundary_list[1]
				nset_end_degree3_dict[Boundary_nset]=Boundary_list[2]
			else:
				print("error in boundary input")
	elif keycase==24:
		pass
	i=i+1
f.close()
g=open("D:\\output4.txt","w")
g.write("Bridge\n")
sum_node_id=0
sum_element_id_T3D2=0
sum_element_id_S4R=0
sum_element_id_B31=0
sum_element_id_C3D8R=0
for tie_surface_list in tie_nset_list:
	iy=0
	slave_surface=tie_surface_list[0]
	slave_surface_nset=surface_nset_dict[slave_surface]
	master_surface=tie_surface_list[1]
	master_surface_nset=surface_nset_dict[master_surface]
	print(len(nset_node_dict[slave_surface_nset]))
	print(len(nset_instance_dict[master_surface_nset]))
	print(slave_surface)
	print(master_surface)
	if len(nset_node_dict[slave_surface_nset])==1:
		help_iy=len(nset_instance_dict[master_surface_nset])
	elif len(nset_instance_dict[master_surface_nset])==1:
		help_iy=len(nset_node_dict[slave_surface_nset])
	else:
		help_iy=len(nset_node_dict[slave_surface_nset])
	while iy<len(nset_node_dict[slave_surface_nset]):
		if len(nset_node_dict[slave_surface_nset])==1:
			slave_surface_node_list=nset_node_dict[slave_surface_nset][0]
			slave_surface_instance=nset_instance_dict[slave_surface_nset][0]
		else:
			slave_surface_node_list=nset_node_dict[slave_surface_nset][iy]
			slave_surface_instance=nset_instance_dict[slave_surface_nset][iy]
		if len(nset_instance_dict[master_surface_nset])==1:
			master_surface_node_list=nset_node_dict[master_surface_nset][0]
			master_surface_instance=nset_instance_dict[master_surface_nset][0]
		else:
			master_surface_node_list=nset_node_dict[master_surface_nset][iy]
			master_surface_instance=nset_instance_dict[master_surface_nset][iy]
		slave_length=len(slave_surface_node_list)
		master_length=len(master_surface_node_list)
		if slave_length!=master_length:
			print("error:the number of the tied points is different")
		else:
			iteration1=0
			while iteration1<slave_length:
				if master_surface_nset=="m_Set-92":
					if slave_surface_nset=="s_Set-92":
						if iteration1==0:
							slave_id=19
							master_id=7
						elif iteration1==1:
							slave_id=20
							master_id=6
						else:
							slave_id=slave_surface_node_list[iteration1]
							master_id=master_surface_node_list[iteration1]
					else:
						slave_id=slave_surface_node_list[iteration1]
						master_id=master_surface_node_list[iteration1]
				else:
					slave_id=slave_surface_node_list[iteration1]
					master_id=master_surface_node_list[iteration1]
				if round(instance_node_dict[slave_surface_instance][slave_id-1][1],7) != round(instance_node_dict[master_surface_instance][master_id-1][1],7):
					print("error:the points location is different")
					print(instance_node_dict[slave_surface_instance][slave_id-1][1])
					print(instance_node_dict[master_surface_instance][master_id-1][1])
					print(slave_id)
					print(master_id)
					print(iteration1)
					print(slave_surface_instance)
					print(master_surface_instance)
					print(master_surface_nset)
					print(slave_surface_nset)
				elif round(instance_node_dict[slave_surface_instance][slave_id-1][2],7) != round(instance_node_dict[master_surface_instance][master_id-1][2],7):
					print("error:the points location is different")
					print(instance_node_dict[slave_surface_instance][slave_id-1][2])
					print(instance_node_dict[master_surface_instance][master_id-1][2])
					print(slave_id)
					print(master_id)
					print(iteration1)
					print(slave_surface_instance)
					print(master_surface_instance)
					print(master_surface_nset)
					print(slave_surface_nset)
				elif round(instance_node_dict[slave_surface_instance][slave_id-1][3],7) != round(instance_node_dict[master_surface_instance][master_id-1][3]):
					print("error:the points location is different")
					print(instance_node_dict[slave_surface_instance][slave_id-1][3])
					print(instance_node_dict[master_surface_instance][master_id-1][3])
					print(slave_id)
					print(master_id)
					print(iteration1)
					print(slave_surface_instance)
					print(master_surface_instance)
				else:
					pass
				instance_node_dict[slave_surface_instance][slave_id-1][4].append([master_surface_instance,master_id])
				instance_node_dict[master_surface_instance][master_id-1][4].append([slave_surface_instance,slave_id])
				iteration1=iteration1+1
		iy=iy+1
help1=0
for every_instance in instance_list:
	this_instance_node_list_old=instance_node_dict[every_instance]
	this_instance_node_list_new=[]
	for this_instance_node in this_instance_node_list_old:
		if len(this_instance_node[4])>0:
			for tie_part in this_instance_node[4]:
				if len(instance_node_dict[tie_part[0]][tie_part[1]-1])==6:
					instance_node_dict[tie_part[0]][tie_part[1]-1][5]=1
					instance_node_dict[tie_part[0]][tie_part[1]-1][0]=this_instance_node[0]+sum_node_id
			if this_instance_node[5]==0:
				new_node_id=this_instance_node[0]+sum_node_id
				new_node_x=this_instance_node[1]
				new_node_y=this_instance_node[2]
				new_node_z=this_instance_node[3]
				this_instance_node_list_new.append([new_node_id,new_node_x,new_node_y,new_node_z])
			elif this_instance_node[5]==1:
				sum_node_id=sum_node_id-1
				help1=help1+1
				new_node_id=this_instance_node[0]
				new_node_x=this_instance_node[1]
				new_node_y=this_instance_node[2]
				new_node_z=this_instance_node[3]
				this_instance_node_list_new.append([new_node_id,new_node_x,new_node_y,new_node_z])
		else:
			new_node_id=this_instance_node[0]+sum_node_id
			new_node_x=this_instance_node[1]
			new_node_y=this_instance_node[2]
			new_node_z=this_instance_node[3]
			this_instance_node_list_new.append([new_node_id,new_node_x,new_node_y,new_node_z])
	instance_node_dict[every_instance]=this_instance_node_list_new
	if every_instance=="Part-RiverBank-2":
		print(this_instance_node_list_new)
	part_name_1=instance_part_dict[every_instance]
	ele_type_1=part_eletype_dict[part_name_1]
	part_element_list_1=part_element_dict[part_name_1]
	instance_element_list_1=[]
	if ele_type_1=="T3D2":
		for part_element_1 in part_element_list_1:
			new_element_id=int(part_element_1[0])+sum_element_id_T3D2
			new_element_data_list=[new_element_id]
			for new_element_node_id in part_element_1[1:]:
				new_element_data_list.append(this_instance_node_list_new[int(new_element_node_id)-1][0])
			instance_element_list_1.append(new_element_data_list)
		sum_element_id_T3D2=sum_element_id_T3D2+len(part_element_list_1)
	elif ele_type_1=="S4R":
		for part_element_1 in part_element_list_1:
			new_element_id=int(part_element_1[0])+sum_element_id_S4R
			new_element_data_list=[new_element_id]
			for new_element_node_id in part_element_1[1:]:
				new_element_data_list.append(this_instance_node_list_new[int(new_element_node_id)-1][0])
			instance_element_list_1.append(new_element_data_list)
		sum_element_id_S4R=sum_element_id_S4R+len(part_element_list_1)
	elif ele_type_1=="C3D8R":
		for part_element_1 in part_element_list_1:
			new_element_id=int(part_element_1[0])+sum_element_id_C3D8R
			new_element_data_list=[new_element_id]
			for new_element_node_id in part_element_1[1:]:
				new_element_data_list.append(this_instance_node_list_new[int(new_element_node_id)-1][0])
			instance_element_list_1.append(new_element_data_list)
		sum_element_id_C3D8R=sum_element_id_C3D8R+len(part_element_list_1)
	elif ele_type_1=="B31":
		for part_element_1 in part_element_list_1:
			new_element_id=int(part_element_1[0])+sum_element_id_B31
			new_element_data_list=[new_element_id]
			for new_element_node_id in part_element_1[1:]:
				new_element_data_list.append(this_instance_node_list_new[int(new_element_node_id)-1][0])
			instance_element_list_1.append(new_element_data_list)
		sum_element_id_B31=sum_element_id_B31+len(part_element_list_1)
	else:
		print("error:unknown element type")
	instance_element_dict[every_instance]=instance_element_list_1
	sum_node_id=sum_node_id+len(this_instance_node_list_old)
print(instance_element_dict["Part-Cable250-2"])
NUMNP=sum_node_id
NUMEG=4
NLCASE=0
MODEX=1
g.write(str(NUMNP)+" "+str(NUMEG)+" "+str(NLCASE)+" "+str(MODEX)+"\n")
id0=1
boundary_id_list=[]
start_degree1_list=[]
end_degree1_list=[]
start_degree2_list=[]
end_degree2_list=[]
start_degree3_list=[]
end_degree3_list=[]
for data_3 in Boundary_nset_list:
	nset_3=data_3
	print(nset_3)
	help_hehe=0
	for instance_4 in nset_instance_dict[nset_3]:
		print(instance_4)
		start_degree=nset_start_degree1_dict[nset_3]
		end_degree=nset_end_degree1_dict[nset_3]
		for data_4 in nset_node_dict[nset_3][help_hehe]:
			#print(data_4)
			#print(instance_node_dict[instance_4])
			data_11=instance_node_dict[instance_4][data_4-1]
			boundary_id_list.append(data_11[0])
			start_degree1_list.append(start_degree)
			end_degree1_list.append(end_degree)
		help_hehe=help_hehe+1
for data_3 in Boundary_nset_list:
	nset_3=data_3
	help_hehe=0
	for instance_4 in nset_instance_dict[nset_3]:
		start_degree=nset_start_degree2_dict[nset_3]
		end_degree=nset_end_degree2_dict[nset_3]
		for data_4 in nset_node_dict[nset_3][help_hehe]:
			start_degree2_list.append(start_degree)
			end_degree2_list.append(end_degree)
		help_hehe=help_hehe+1
for data_3 in Boundary_nset_list:
	nset_3=data_3
	help_hehe=0
	for instance_4 in nset_instance_dict[nset_3]:
		start_degree=nset_start_degree3_dict[nset_3]
		end_degree=nset_end_degree3_dict[nset_3]
		for data_4 in nset_node_dict[nset_3][help_hehe]:
			start_degree3_list.append(start_degree)
			end_degree3_list.append(end_degree)
		help_hehe=help_hehe+1
node_output_list=[]
for instance_3 in instance_list:
	for node_list_3 in instance_node_dict[instance_3]:
		if node_list_3[0]<id0:
			pass
		else:
			id0=node_list_3[0]
			N=id0
			ID1=0
			ID2=0
			ID3=0
			X=node_list_3[1]
			Y=node_list_3[2]
			Z=node_list_3[3]
			node_output_list.append([N,ID1,ID2,ID3,X,Y,Z])
num_degree=0
for boundary_id in boundary_id_list:
	start_degree=int(start_degree1_list[num_degree])
	end_degree=int(end_degree1_list[num_degree])
	num_degree=num_degree+1
	iteration2=start_degree
	while iteration2<end_degree+1:
		node_output_list[boundary_id-1][iteration2]=1
		iteration2=iteration2+1
num_degree=0
for boundary_id in boundary_id_list:
	start_degree=int(start_degree2_list[num_degree])
	end_degree=int(end_degree2_list[num_degree])
	num_degree=num_degree+1
	iteration2=start_degree
	while iteration2<end_degree+1:
		node_output_list[boundary_id-1][iteration2]=1
		iteration2=iteration2+1
num_degree=0
for boundary_id in boundary_id_list:
	start_degree=int(start_degree3_list[num_degree])
	end_degree=int(end_degree3_list[num_degree])
	num_degree=num_degree+1
	iteration2=start_degree
	while iteration2<end_degree+1:
		node_output_list[boundary_id-1][iteration2]=1
		iteration2=iteration2+1
for output_node_line in node_output_list:
	N=output_node_line[0]
	ID1=output_node_line[1]
	ID2=output_node_line[2]
	ID3=output_node_line[3]
	X=output_node_line[4]
	Y=output_node_line[5]
	Z=output_node_line[6]
	g.write(str(N)+" "+str(ID1)+" "+str(ID2)+" "+str(ID3)+" "+str(X)+" "+str(Y)+" "+str(Z)+"\n")
NPAR1=1
NPAR2=sum_element_id_T3D2
NPAR3=0
T3D2_output_data=[]
MN_MTYP_dict={"test":0}
for data_6 in instance_list:
	part_name=instance_part_dict[data_6]
	material_name1=part_section_material_dict[part_name]
	ele_type_2=part_eletype_dict[part_name]
	if ele_type_2=="T3D2":
		for data_7 in instance_element_dict[data_6]:
			M=data_7[0]
			II=data_7[1]
			JJ=data_7[2]
			MN=material_name1
			AREA=float(part_material_para_dict[part_name][0])
			T3D2_output_data.append([M,II,JJ,MN,AREA])
	else:
		pass
help_list1=[]
for help_x in T3D2_output_data:
	help_list1.append(help_x[3])
help_list2=[]
for help_y in T3D2_output_data:
	help_list2.append(help_y[4])
for mat1 in set(help_list1):
	for mat2 in set(help_list2):
		NPAR3=NPAR3+1
		MN_MTYP_dict[(mat1,mat2)]=NPAR3
g.write(str(NPAR1)+" "+str(NPAR2)+" "+str(NPAR3)+"\n")
for mtyp1 in MN_MTYP_dict.keys():
	if mtyp1 =="test":
		pass
	else:
		N=MN_MTYP_dict[mtyp1]
		E=float(material_E_dict[mtyp1[0]])
		poisson=float(material_poisson_dict[mtyp1[0]])
		density=float(material_density_dict[mtyp1[0]])
		AREA=float(mtyp1[1])
		g.write(str(N)+" "+str(E)+" "+str(poisson)+" "+str(AREA)+" "+str(density)+"\n")
for data_8 in T3D2_output_data:
	M=data_8[0]
	II=data_8[1]
	JJ=data_8[2]
	MN=data_8[3:5]
	MTYP=MN_MTYP_dict[MN[0],MN[1]]
	g.write(str(M)+" "+str(II)+" "+str(JJ)+" "+str(MTYP)+"\n")
NPAR1=2
NPAR2=sum_element_id_C3D8R
NPAR3=0
T3D2_output_data=[]
MN_MTYP_dict={"test":0}
for data_6 in instance_list:
	part_name=instance_part_dict[data_6]
	material_name1=part_section_material_dict[part_name]
	ele_type_2=part_eletype_dict[part_name]
	if ele_type_2=="C3D8R":
		for data_7 in instance_element_dict[data_6]:
			M=data_7[0]
			II=data_7[1]
			JJ=data_7[2]
			MN=material_name1
			T3D2_output_data.append([M,II,JJ,MN])
	else:
		pass
help_list1=[]
for help_x in T3D2_output_data:
	help_list1.append(help_x[3])
for mat1 in set(help_list1):
	NPAR3=NPAR3+1
	MN_MTYP_dict[mat1]=NPAR3
g.write(str(NPAR1)+" "+str(NPAR2)+" "+str(NPAR3)+"\n")
for mtyp1 in MN_MTYP_dict.keys():
	if mtyp1 =="test":
		pass
	else:
		N=MN_MTYP_dict[mtyp1]
		E=float(material_E_dict[mtyp1])
		poisson=float(material_poisson_dict[mtyp1])
		denstity=float(material_poisson_dict[mtyp1])
		g.write(str(N)+" "+str(E)+" "+str(poisson)+" "+str(density)+"\n")
for data_8 in T3D2_output_data:
	M=data_8[0]
	II=data_8[1]
	JJ=data_8[2]
	MN=data_8[3]
	MTYP=MN_MTYP_dict[MN]
	g.write(str(M)+" "+str(II)+" "+str(JJ)+" "+str(MTYP)+"\n")
NPAR1=2
NPAR2=sum_element_id_S4R
NPAR3=0
T3D2_output_data=[]
MN_MTYP_dict={"test":0}
for data_6 in instance_list:
	part_name=instance_part_dict[data_6]
	material_name1=part_section_material_dict[part_name]
	ele_type_2=part_eletype_dict[part_name]
	if ele_type_2=="S4R":
		for data_7 in instance_element_dict[data_6]:
			M=data_7[0]
			II=data_7[1]
			JJ=data_7[2]
			MN=material_name1
			WIDTH=float(part_material_para_dict[part_name][0])
			T3D2_output_data.append([M,II,JJ,MN,WIDTH])
	else:
		pass
help_list1=[]
for help_x in T3D2_output_data:
	help_list1.append(help_x[3])
help_list2=[]
for help_y in T3D2_output_data:
	help_list2.append(help_y[4])
for mat1 in set(help_list1):
	for mat2 in set(help_list2):
		NPAR3=NPAR3+1
		MN_MTYP_dict[(mat1,mat2)]=NPAR3
g.write(str(NPAR1)+" "+str(NPAR2)+" "+str(NPAR3)+"\n")
for mtyp1 in MN_MTYP_dict.keys():
	if mtyp1 =="test":
		pass
	else:
		N=MN_MTYP_dict[mtyp1]
		E=float(material_E_dict[mtyp1[0]])
		poisson=float(material_poisson_dict[mtyp1[0]])
		density=float(material_density_dict[mtyp1[0]])
		WIDTH=float(mtyp1[1])
		g.write(str(N)+" "+str(E)+" "+str(poisson)+" "+str(density)+" "+str(WIDTH)+"\n")
for data_8 in T3D2_output_data:
	M=data_8[0]
	II=data_8[1]
	JJ=data_8[2]
	MN=data_8[3:5]
	MTYP=MN_MTYP_dict[MN[0],MN[1]]
	g.write(str(M)+" "+str(II)+" "+str(JJ)+" "+str(MTYP)+"\n")
NPAR1=2
NPAR2=sum_element_id_B31
NPAR3=0
T3D2_output_data=[]
MN_MTYP_dict={"test":0}
for data_6 in instance_list:
	part_name=instance_part_dict[data_6]
	material_name1=part_section_material_dict[part_name]
	ele_type_2=part_eletype_dict[part_name]
	if ele_type_2=="B31":
		for data_7 in instance_element_dict[data_6]:
			M=data_7[0]
			II=data_7[1]
			JJ=data_7[2]
			MN=material_name1
			data11=float(part_material_para_dict[part_name][0][0])
			data12=float(part_material_para_dict[part_name][0][1])
			data13=float(part_material_para_dict[part_name][0][2])
			data14=float(part_material_para_dict[part_name][0][4])
			data21=float(part_material_para_dict[part_name][1][0])
			data22=float(part_material_para_dict[part_name][1][1])
			data23=float(part_material_para_dict[part_name][1][2])
			T3D2_output_data.append([M,II,JJ,MN,data11,data12,data13,data14,data21,data22,data23])
	else:
		pass
help_list1=[]
for help_x in T3D2_output_data:
	help_list1.append(help_x[3])
help_list2=[]
for help_y in T3D2_output_data:
	help_list2.append(help_y[4])
help_2=[]
for help_y2 in T3D2_output_data:
	help_2.append(help_y2[5])
help_3=[]
for help_y2 in T3D2_output_data:
	help_3.append(help_y2[6])
help_4=[]
for help_y2 in T3D2_output_data:
	help_4.append(help_y2[7])
help_5=[]
for help_y2 in T3D2_output_data:
	help_5.append(help_y2[8])
help_6=[]
for help_y2 in T3D2_output_data:
	help_6.append(help_y2[9])
help_7=[]
for help_y2 in T3D2_output_data:
	help_7.append(help_y2[10])
for mat1 in set(help_list1):
	for mat2 in set(help_list2):
		for mat3 in set(help_2):
			for mat4 in set(help_3):
				for mat5 in set(help_4):
					for mat6 in set(help_5):
						for mat7 in set(help_6):
							for mat8 in set(help_7):
								NPAR3=NPAR3+1
								MN_MTYP_dict[(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)]=NPAR3
g.write(str(NPAR1)+" "+str(NPAR2)+" "+str(NPAR3)+"\n")
for mtyp1 in MN_MTYP_dict.keys():
	if mtyp1 =="test":
		pass
	else:
		N=MN_MTYP_dict[mtyp1]
		E=float(material_E_dict[mtyp1[0]])
		poisson=float(material_poisson_dict[mtyp1[0]])
		data11=float(mtyp1[1])
		data12=float(mtyp1[2])
		data13=float(mtyp1[3])
		data14=float(mtyp1[4])
		data21=float(mtyp1[5])
		data22=float(mtyp1[6])
		data23=float(mtyp1[7])
		g.write(str(N)+" "+str(E)+" "+str(poisson)+" "+str(data11)+" "+str(data12)+" "+str(data13)+" "+str(data14)+" "+str(data21)+" "+str(data22)+" "+str(data23)+" "+"\n")
for data_8 in T3D2_output_data:
	M=data_8[0]
	II=data_8[1]
	JJ=data_8[2]
	MN=data_8[3:11]
	MTYP=MN_MTYP_dict[MN[0],MN[1],MN[2],MN[3],MN[4],MN[5],MN[6],MN[7]]
	g.write(str(M)+" "+str(II)+" "+str(JJ)+" "+str(MTYP)+"\n")
g.close()