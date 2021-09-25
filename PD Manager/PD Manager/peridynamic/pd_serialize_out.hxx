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
				}
				void		OutputHwAscii(double curtime)
				{
					const TPdModel& pdModel = *m_p_PdModel;
					int width = 20;

					fprintf(fout, "$TIME\t = %15.10f sec\n", curtime);

					fprintf(fout, "$BINDING\t = NODE\n");
					fprintf(fout, "$COLUMN_INFO\t = ENTITY_ID\n");
					fprintf(fout, "$RESULT_TYPE\t = Displacement(v), Rotation(v)\n");

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

						fprintf(fout, "%20.10E", node.Displacement().rx());
						fprintf(fout, "%20.10E", node.Displacement().ry());
						fprintf(fout, "%20.10E", node.Displacement().rz());
						
						fprintf(fout, "\n");
					}	

					/************************************************************************/
					/* 按照有限元单元输出物理量                                             */
					/************************************************************************/
					fprintf(fout, "$BINDING\t = ELEMENT\n");
					fprintf(fout, "$COLUMN_INFO\t = ENTITY_ID\n");
					fprintf(fout, "$RESULT_TYPE\t = DamageIndex(s)\n");

					const set<int> femElemIds = pdModel.PdMeshCore().GetElementIdsByAll();
					for (int eid : femElemIds)
					{
						const TPdElement& element = pdModel.PdMeshCore().Element(eid);
						fprintf(fout, "%15d", element.Id() + 1);

						fprintf(fout, "%20.10E", element.DamageIndex());

						fprintf(fout, "\n");
					}
				}
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif