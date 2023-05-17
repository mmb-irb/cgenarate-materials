#!/usr/bin/env python
# coding: utf-8
# In[3]:


import dnatraj as dnat
from dnatraj import duplex as dx


# In[4]:


import pickle


# In[5]:


import mdtraj as mdt
from mdplus.multiscale import Glimps
import numpy as np
	
# In[1]:


names = ["1naj",]

for i in range(len(names)):


	# In[2]:


	referencetop = 'input/' + names[i] + '.pdb'
	inputtraj = 'input/' + names[i] + 'aligned.mdcrd'
	outputtop = 'output/' + names[i] + 'predicted.pdb'
	outputtraj = 'output/' + names[i] + 'predicted.mdcrd'

	print(names[i])



	# In[6]:


	# bp=40
	backbone_nmers=10

	#bases = bp*2
	#firsttestbp = bp-backbone_nmers-1
	#trainnmers = firsttestbp-backbone_nmers

	#trainbp = (firsttestbp-1)
	#trainbases = 2*trainbp

	#testbp = backbone_nmers
	#testbases = 2*testbp


	# # Rebuild full oligo

	# In[7]:


	with open("transformer_backbone.pickle", "rb") as input_file:
	    transformer_backbone=pickle.load(input_file)
	with open("transformerA.pickle", "rb") as input_file:
	    transformerA=pickle.load(input_file)
	with open("transformerC.pickle", "rb") as input_file:
	    transformerC=pickle.load(input_file)


	# In[ ]:





	# In[8]:



	mypdb = mdt.load(referencetop)


	# In[9]:


	sel = mypdb.topology.select('not type H')
	mynoHpdb = mypdb.atom_slice(sel)

	mynoHpdb.save(outputtop)


	# In[10]:


	sel = mypdb.topology.select('name =~ "C1."')
	myc1pdb = mypdb.atom_slice(sel)


	# In[11]:


	cg_test_backbone = mdt.load(inputtraj, top=myc1pdb.top)[::10]
	dx_test_backbone = dx.ComplementaryDuplex(cg_test_backbone)
	lendx=len(dx_test_backbone)


	# In[12]:


	cg_test_backbone


	# In[13]:


	dx_test_backbone.sequence#,dx_test_backbone.sequence==Sequences[simnum-1]


	# In[14]:


	# cg_test_backbone.save("AtomAAAA_skip.mdcrd")


	# In[15]:


	#Topology with extra phosphates for full reconstruction
	auxtop=mynoHpdb.top.copy()

	for name in ("OP2", "OP1", "P"):
	    auxtop.insert_atom(name=name,
		               element=next(auxtop.atoms_by_name(name)).element,
		               residue=auxtop._residues[0],
		               index=0,
		               rindex=0)
	    
	for name in ("OP2", "OP1", "P"):
	    auxtop.insert_atom(name=name,
		               element=next(auxtop.atoms_by_name(name)).element,
		               residue=auxtop._residues[lendx],
		               index=auxtop._residues[lendx]._atoms[0].index,
		               rindex=0) 

	auxpdb=mdt.Trajectory(np.zeros((1,auxtop._numAtoms,3)),auxtop)


	# ## Rebuild backbone

	# ### Rebuild10mers

	# In[16]:



	dx_nmer_pred_backbone = []

	sel = auxpdb.topology.select('name =~ ".*\'.*|.*P.*"') 
	dx_mypdb_backbone = dx.ComplementaryDuplex(auxpdb.atom_slice(sel))

	for i in range(lendx // backbone_nmers):

	    dx_backbone = dx_test_backbone[backbone_nmers * i:backbone_nmers * (i + 1)]
	    
	    xfg_predtest_backbone = transformer_backbone.transform(
		dx_backbone.traj.xyz)

	    dx_mypdbs = dx_mypdb_backbone[backbone_nmers * i:backbone_nmers * (i + 1)]
	    dx_pred = mdt.Trajectory(xfg_predtest_backbone, dx_mypdbs.traj.topology)
	    dx_nmer_pred_backbone.append(dx.ComplementaryDuplex(dx_pred))
	    print(i)

	if (lendx % backbone_nmers != 0):

	    dx_backbone = dx_test_backbone[lendx - backbone_nmers:lendx]

	    xfg_predtest_backbone = transformer_backbone.transform(
		dx_backbone.traj.xyz)

	    dx_mypdbs = dx_mypdb_backbone[lendx - backbone_nmers:lendx]
	    dx_pred = mdt.Trajectory(xfg_predtest_backbone, dx_mypdbs.traj.topology)
	    dx_nmer_pred_backbone.append(dx.ComplementaryDuplex(dx_pred))
	    print("End")


	# ### Join 10mers

	# In[17]:



	cut = backbone_nmers - lendx % backbone_nmers
	dx_last_pred_backbone = dx_nmer_pred_backbone[-1][cut:backbone_nmers]
	dx_pred_backbone = dx.sstack(dx_nmer_pred_backbone[:-1] +
		                    [dx_last_pred_backbone])
	fg_predtest_backbone = dx_pred_backbone.traj


	# In[18]:


	# fg_predtest_backbone.save_pdb('test_backbone.pdb')


	# ## Rebuild Adenine from backbone 

	# Rebuilt Adenine in test from "fg_predtest_backbone"

	# In[19]:




	Aindexes_test_backbone = [
	   res.index for res in cg_test_backbone.topology._residues
	   if res.name[:2] == 'DA'
	]

	cg_test_Afrombackbone = [None] * len(Aindexes_test_backbone)

	for i, index in enumerate(Aindexes_test_backbone):
	   vmdsel = 'resid ' + str(index)
	   sel = fg_predtest_backbone.topology.select(vmdsel)
	   tempA = fg_predtest_backbone.atom_slice(sel)

	   vmdsel = 'resid ' + str(lendx * 2 - 1 - index)
	   sel = fg_predtest_backbone.topology.select(vmdsel)
	   tempT = fg_predtest_backbone.atom_slice(sel)

	   cg_test_Afrombackbone[i] = tempA.stack(tempT, True)


	# In[20]:




	testtop=auxtop

	fg_predtest_Afrombackbone = [None]*len(Aindexes_test_backbone)

	for i, index in enumerate(Aindexes_test_backbone):
	   
	   vmdsel = 'resid ' + str(index)
	   sel = testtop.select(vmdsel)
	   topA = testtop.subset(sel)

	   vmdsel = 'resid ' + str(lendx*2-1-index)
	   sel = testtop.select(vmdsel)
	   topT = testtop.subset(sel)
	   
	   topAT = topA.join(topT)
	   
	   xfg_predtest_Afrombackbone = transformerA.transform(cg_test_Afrombackbone[i].xyz)
	   #merge with topology
	   fg_predtest_Afrombackbone[i] = mdt.Trajectory(xfg_predtest_Afrombackbone, topAT) 
	
	print("AT")
	   
	#     sel1 = fg_predtest_Afrombackbone[i].topology.select('name =~ ".*\'.*|.*P.*"')
	#     sel2 = cg_test_Afrombackbone[i].topology.select('all')
	#     for j in range(len(fg_predtest_Afrombackbone[i])):     
	#         fg_predtest_Afrombackbone[i][j].superpose(cg_test_Afrombackbone[i][j],atom_indices=sel1,ref_atom_indices=sel2)
	#    print(index)


	# In[21]:


	# def mdtstack (trajectories, keep_resSeq=True):
	#     #trajectories to stack (same time)
	    
	#     if len(trajectories) < 2:
	#         print('why')
	#         return
	    
	#     finaltraj=trajectories[0].stack(trajectories[1], keep_resSeq=True)
	    
	#     for i in trajectories[2:]:
	#         finaltraj=finaltraj.stack(i, keep_resSeq=True)
		
	#     return finaltraj


	# In[22]:


	# mdtstack(fg_predtest_Afrombackbone, keep_resSeq=True).save_pdb('predtest_Afrombackbone.pdb')


	# ## Rebuild Cytosine from backbone 

	# Rebuilt Cytosine in test from "fg_predtest_backbone"

	# In[23]:




	#select based on original traj
	Cindexes_test_backbone = [res.index for res in cg_test_backbone.topology._residues if res.name[:2] == 'DC'] 

	cg_test_Cfrombackbone = [None]*len(Cindexes_test_backbone)

	for i, index in enumerate(Cindexes_test_backbone):
	   vmdsel = 'resid ' + str(index)
	   sel = fg_predtest_backbone.topology.select(vmdsel)
	   tempC = fg_predtest_backbone.atom_slice(sel)
	   
	   vmdsel = 'resid ' + str(lendx*2-1-index)
	   sel = fg_predtest_backbone.topology.select(vmdsel)
	   tempG = fg_predtest_backbone.atom_slice(sel)
	   
	   cg_test_Cfrombackbone[i]=tempC.stack(tempG,True)


	# In[24]:




	testtop=auxtop

	fg_predtest_Cfrombackbone = [None]*len(Cindexes_test_backbone)

	for i, index in enumerate(Cindexes_test_backbone):
	   
	   vmdsel = 'resid ' + str(index)
	   sel = testtop.select(vmdsel)
	   topC = testtop.subset(sel)

	   vmdsel = 'resid ' + str(lendx*2-1-index)
	   sel = testtop.select(vmdsel)
	   topG = testtop.subset(sel)
	   
	   topCG = topC.join(topG)
	   
	   xfg_predtest_Cfrombackbone = transformerC.transform(cg_test_Cfrombackbone[i].xyz)
	   fg_predtest_Cfrombackbone[i] = mdt.Trajectory(xfg_predtest_Cfrombackbone, topCG)
	   
	print("CG")
	#     sel1 = fg_predtest_Cfrombackbone[i].topology.select('name =~ ".*\'.*|.*P.*"')
	#     sel2 = cg_test_Cfrombackbone[i].topology.select('all')
	#     for j in range(len(fg_predtest_Cfrombackbone[i])):
	#         fg_predtest_Cfrombackbone[i][j].superpose(cg_test_Cfrombackbone[i][j],atom_indices=sel1,ref_atom_indices=sel2)
	#    print(index)


	# In[25]:


	# def mdtstack (trajectories, keep_resSeq=True):
	#     #trajectories to stack (same time)
	    
	#     if len(trajectories) < 2:
	#         print('why')
	#         return
	    
	#     finaltraj=trajectories[0].stack(trajectories[1], keep_resSeq=True)
	    
	#     for i in trajectories[2:]:
	#         finaltraj=finaltraj.stack(i, keep_resSeq=True)
		
	#     return finaltraj


	# In[26]:


	# mdtstack(fg_predtest_Cfrombackbone, keep_resSeq=True).save_pdb('predtest_Cfrombackbone.pdb')


	# ## Rejoin A and C predictions in backbone

	# In[27]:



	tostack=[]
	Astack=iter(fg_predtest_Afrombackbone)
	Tstack=iter(fg_predtest_Afrombackbone[::-1])
	Cstack=iter(fg_predtest_Cfrombackbone)
	Gstack=iter(fg_predtest_Cfrombackbone[::-1])

	for i in range(lendx):
	   if i in Aindexes_test_backbone:
	       tostack.append(dx.ComplementaryDuplex(next(Astack)))
	   elif lendx*2-1-i in Aindexes_test_backbone:
	       tostack.append(dx.ComplementaryDuplex(next(Tstack)).invert())
	   elif i in Cindexes_test_backbone:
	       tostack.append(dx.ComplementaryDuplex(next(Cstack)))
	   elif lendx*2-1-i in Cindexes_test_backbone:
	       tostack.append(dx.ComplementaryDuplex(next(Gstack)).invert())
	   else:
	       print("Error",i)


	# In[28]:


	fg_predtest_0 = dx.sstack(tostack) 


	# In[29]:


	# fg_predtest_0.traj.save_pdb('test_predicted0.pdb')


	# In[30]:


	sel = fg_predtest_0.traj.topology.select(
	    'not ((resid 0 %d) and (name =~ ".*P.*"))' % (lendx))
	fg_predtest = fg_predtest_0.traj.atom_slice(sel)


	# In[31]:


	#  
	# import time
	# print(time.asctime(time.localtime()))


	# In[32]:


	print("Finished")

	fg_predtest.save(outputtraj)
