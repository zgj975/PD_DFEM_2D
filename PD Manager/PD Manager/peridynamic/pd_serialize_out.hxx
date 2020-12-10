/*==============================================================================

Copyright 2017 Dalian University of Technology .
All rights reserved

================================================================================
-- Please append file description informations here --
Write PD results
SerializeHwAscii	:	Output results for HyperView (plot)
SerializeNodesInfo	:	Output results for Origin (curves) 
================================================================================
Date            Name                    Description of Change
2018/04/05		Zheng Guojun			Create
2019/07/31		Zheng Guojun			Modified to TEMPLATE type
$HISTORY$
================================================================================
*/

#ifndef DLUT_SAE_PERIDYNAMIC_PD_SERIALIZE_OUT_HXX_20171224
#define DLUT_SAE_PERIDYNAMIC_PD_SERIALIZE_OUT_HXX_20171224

#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <iomanip>
using namespace std;

#include "pd_database.hxx"

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			class SerializeBase
			{
			public:
				SerializeBase()
				{
					m_p_PdModel = NULL;
					fout = NULL;
				}
				~SerializeBase()
				{
					DetachFile();
				}
			public:
				void		AttachModel(TPdModel& pdModel)
				{
					m_p_PdModel = &pdModel;
				}
				void		AttachFile(string filename)
				{
					fopen_s(&fout, filename.c_str(), "wt");
				}
				void		DetachFile()
				{
					m_p_PdModel = NULL;
					if (fout != NULL)
					{
						fclose(fout);
					}
				}
			protected:
				FILE*		fout;
				TPdModel*	m_p_PdModel;
			};

			class SerializeHwAscii : public SerializeBase
			{
			public:
				void		AttachFile(string filename)
				{
					SerializeBase::AttachFile(filename);

					fprintf(fout, "ALTAIR ASCII FILE\n");
					fprintf(fout, "$TITLE\t = Transient analysis\n");
					fprintf(fout, "$SUBCASE\t = 1\tSubcase\t1\n");
					fprintf(fout, "$BINDING\t = NODE\n"); 
				//	fprintf(fout, "$BINDING\t = ELEMENT\n");
					fprintf(fout, "$COLUMN_INFO\t = ENTITY_ID\n");
					fprintf(fout, "$RESULT_TYPE\t = Displacement(v), DamageIndex(s)\n");
				//	fprintf(fout, "$RESULT_TYPE\t = DamageIndex(s)\n");
				}
				void		OutputHwAscii(double curtime)
				{
					const TPdModel& pdModel = *m_p_PdModel;
					int width = 20;

					fprintf(fout, "$TIME\t = %15.10f sec\n", curtime);

					/************************************************************************/
					/* 按照有限元节点输出物理量                                             */
					/************************************************************************/
					set<int> femNodeIds = pdModel.PdMeshCore().GetNodeIdsByAll();

					for (int nid : femNodeIds)
					{
						const TPdNode& node = pdModel.PdMeshCore().Node(nid);
												
						fprintf(fout, "%15d", node.Id()+1);

						fprintf(fout, "%20.10E", node.Displacement().x());
						fprintf(fout, "%20.10E", node.Displacement().y());
						fprintf(fout, "%20.10E", node.Displacement().z());

						set<int> adjElems = node.AdjElementIds();
						double damageIndex = 0;
						for (int ei : adjElems)
						{
							const TPdElement& element = pdModel.PdMeshCore().Element(ei);
							damageIndex += element.DamageIndex();
						}
						damageIndex /= (int)(adjElems.size());
						fprintf(fout, "%20.10E", damageIndex);

						fprintf(fout, "\n");
					}	

					/*set<int> femEleIds = pdModel.PdMeshCore().GetElementIdsByAll();
					for (int eleId : femEleIds)
					{
						const TPdElement& element = pdModel.PdMeshCore().Element(eleId);
						fprintf(fout, "%15d", element.Id() + 1);
						fprintf(fout, "%20.10E", element.DamageIndex());

						fprintf(fout, "\n");
					}*/
				}
			};

			class SerializeElementsInfo : public SerializeBase
			{
			public:
				void		AttachFile(string filename)
				{
					if (m_p_PdModel->PdMeshCore().GetOutputElementIds().size() > 0)
					{
						SerializeBase::AttachFile(filename);

						const TPdModel& pdModel = *m_p_PdModel;
						fprintf(fout, "%10s", "Time");

						const set<int>& opEleIds = pdModel.PdMeshCore().GetOutputElementIds();
						for (int eid : opEleIds)
						{
							const TElementBase& ele = pdModel.PdMeshCore().Element(eid);

							fprintf(fout, "Ele_%15d_DX", ele.Id()+1);
							fprintf(fout, "Ele_%15d_DY", ele.Id()+1);
							fprintf(fout, "Ele_%15d_DZ", ele.Id()+1);
						}
						fprintf(fout, "\n");
					}
				}
				void		OutputElementInfo(double curtime)
				{
					if (fout != NULL)
					{
						const TPdModel& pdModel = *m_p_PdModel;
						int width = 20;
						int precision = width / 2;
						const set<int>& opEleIds = pdModel.PdMeshCore().GetOutputElementIds();

						fprintf(fout, "%15.10f", curtime);
						for (int eid : opEleIds)
						{
							const TPdElement& element = pdModel.PdMeshCore().Element(eid);
							TDisplacement nodeDis = pdModel.PdMeshCore().Node(element.NodeId(0)).Displacement() +
								pdModel.PdMeshCore().Node(element.NodeId(1)).Displacement() +
								pdModel.PdMeshCore().Node(element.NodeId(2)).Displacement() +
								pdModel.PdMeshCore().Node(element.NodeId(3)).Displacement();
							nodeDis /= 4;

							fprintf(fout, "%20.10f", nodeDis.x());						
							fprintf(fout, "%20.10f", element.DisplaceInElement(0, 0).x());

							fprintf(fout, "%20.10f", nodeDis.y());
							fprintf(fout, "%20.10f", element.DisplaceInElement(0, 0).y());

							fprintf(fout, "%20.10f", nodeDis.z());
							fprintf(fout, "%20.10f", element.DisplaceInElement(0, 0).z());
						}
						fprintf(fout, "\n");
					}
				}
			};

			class SerializeNodesInfo : public SerializeBase
			{
			public:
				void		AttachFile(string filename)
				{
					if (m_p_PdModel->PdMeshCore().GetOutputNodeIds().size() > 0)
					{
						SerializeBase::AttachFile(filename);

						fprintf(fout, "%10s\t", "Time");
						const TPdModel& pdModel = *m_p_PdModel;
						const set<int>& opNodeIds = pdModel.PdMeshCore().GetOutputNodeIds();
						for (int nid : opNodeIds)
						{
							const TPdNode& femNode = pdModel.PdMeshCore().Node(nid);

							fprintf(fout, "Node_%d_DX\t\t", femNode.IdGlobal());
							fprintf(fout, "Node_%d_DY\t\t", femNode.IdGlobal());
							fprintf(fout, "Node_%d_DZ\t\t", femNode.IdGlobal());
						}
						fprintf(fout, "\n");
					}

				}
				void		OutputNodesInfo(double curtime)
				{
					if (fout != NULL)
					{
						fprintf(fout, "%15.10f", curtime);
						const TPdModel& pdModel = *m_p_PdModel;
						const set<int>& opNodeIds = pdModel.PdMeshCore().GetOutputNodeIds();
						for (int nid : opNodeIds)
						{
							const TPdNode& femNode = pdModel.PdMeshCore().Node(nid);
							
							fprintf(fout, "%20.10E", femNode.InnerForce().x());
							fprintf(fout, "%20.10E", femNode.InnerForce().y());
							fprintf(fout, "%20.10E", femNode.InnerForce().z());
						}
						fprintf(fout, "\n");
					}
				}
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif