/*All rights reserved

== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
--Please append file description informations here --
Read mesh informations from a LSDYNA file (.key) to PD database
1. Read *NODE
2. Read *ELEMENT 
3. Insert NODES to PD database 
4. Read other information from .key file, such as FORCE/VELOCITY/DISPLACEMENT/INITIAL CRACK etc. 
== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
Date			    Name                    Description of Change
2018 / 04 / 05		Zheng Guojun			Create
2019 / 07 / 31		Zheng Guojun			Modified the GLOBAL <-> PD NODE LOCAL ID
$HISTORY$
== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
*/

#ifndef DLUT_SAE_PERIDYNAMIC_PD_SERIALIZE_IN_LSDYNA_HXX_20171224
#define DLUT_SAE_PERIDYNAMIC_PD_SERIALIZE_IN_LSDYNA_HXX_20171224

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
			class SerializeLsdyna
			{
			public:
				SerializeLsdyna() { /*Do nothing*/ }
			public:				
				void		ReadLsdynaFile(string filename, TPdModel& pdModel)
				{
					double start = clock();

					readMeshData(filename, pdModel);
					readCalculateInfo(filename, pdModel);
					readBoundaryInfo(filename, pdModel);

					double total_time = (clock() - start) / 1000;
					cout << "ReadLsdynaFile(): \t\t" << total_time << endl;
				}
				void		WriteLsdynaFile(string filename, TPdModel& pdModel)
				{
					FILE* fout;
					fopen_s(&fout, filename.c_str(), "wt");
					if (fout != NULL)
					{
						fprintf(fout, "*KEYWORD\n");
						fprintf(fout, "*NODE\n");
						set<int> femNodeIds = pdModel.PdMeshCore().GetNodeIdsByAll();
						for (int nid : femNodeIds)
						{
							const TPdNode& node = pdModel.PdMeshCore().Node(nid);
							fprintf(fout, "%8d%16f%16f%16f\n", node.Id()+1, node.Coordinate().x(), node.Coordinate().y(), node.Coordinate().z());
						}
												
						fprintf(fout, "*ELEMENT_SHELL\n");
						set<int> femElementIds = pdModel.PdMeshCore().GetElementIdsByAll();
						for (int eid : femElementIds)
						{
							const TPdElement& element = pdModel.PdMeshCore().Element(eid);
							if (element.MeshElementType() == TRIANGLE_ELEMENT)
							{
								fprintf(fout, "%8d%8d%8d%8d%8d%8d\n", element.Id() + 1, element.PartId(), element.NodeId(0) + 1, element.NodeId(1) + 1, element.NodeId(2) + 1, element.NodeId(2) + 1);
							}
							else if (element.MeshElementType() == QUADRANGLE_ELEMENT)
							{
								fprintf(fout, "%8d%8d%8d%8d%8d%8d\n", element.Id() + 1, element.PartId(), element.NodeId(0) + 1, element.NodeId(1) + 1, element.NodeId(2) + 1, element.NodeId(3) + 1);
							}
						}

						fprintf(fout, "*END\n");
						fclose(fout);
					}
				}
			private:
				void		readMeshData(string filename, TPdModel& pdModel)
				{
					pdModel.Initialize();
					map_fem_nid_global_to_fem_nid_local.clear();

					ifstream fin(filename, ios::in);
					if (!fin.is_open())
						return;

					string line;
					int data_type = 0;

					while (!fin.eof())
					{
						getline(fin, line);
						if (line.length() < 1)
							continue;
						if (line.substr(0, 1) == string("$"))
							continue;

						if (string("*NODE") == line)
						{
							streampos pos = fin.tellg();
							while (getline(fin, line))
							{
								if (line.substr(0, 1) == string("*") ||
									line.substr(0, 1) == string("$"))
								{
									fin.seekg(pos);
									break;
								}
								int nid = ParseString<int>(line.substr(0, 8));
								double x = ParseString<double>(line.substr(9, 16));
								double y = ParseString<double>(line.substr(25, 16));
								double z = ParseString<double>(line.substr(41, 16));

								pdModel.PdMeshCore().InsertNode(Vector3d(x, y, z), nid);
								map_fem_nid_global_to_fem_nid_local.insert(pair<int, int>(nid, pdModel.PdMeshCore().NodeCount()-1));
								
								pos = fin.tellg();
							}

							continue;
						}
						else if (string("*ELEMENT_SHELL") == line)
						{
							streampos pos = fin.tellg();
							while (getline(fin, line))
							{
								if (line.substr(0, 1) == string("*") ||
									line.substr(0, 1) == string("$"))
								{
									fin.seekg(pos);
									break;
								}

								int eid = ParseString<int>(line.substr(0, 8));
								int pid = ParseString<int>(line.substr(8, 8));
								int n1 = ParseString<int>(line.substr(16, 8));
								int n2 = ParseString<int>(line.substr(24, 8));
								int n3 = ParseString<int>(line.substr(32, 8));
								int n4 = ParseString<int>(line.substr(40, 8));

								vector<int> nids;
								nids.push_back(map_fem_nid_global_to_fem_nid_local[n1]);
								nids.push_back(map_fem_nid_global_to_fem_nid_local[n2]);
								nids.push_back(map_fem_nid_global_to_fem_nid_local[n3]);

								if (n3 != n4)
								{
									nids.push_back(map_fem_nid_global_to_fem_nid_local[n4]);
								}

								pdModel.PdMeshCore().InsertElement(nids, eid, pid);
								map_fem_eid_global_to_fem_eid_local.insert(pair<int, int>(eid, pdModel.PdMeshCore().ElementCount() - 1));

								pos = fin.tellg();
							}

							continue;
						}
					}
					fin.close();
				}
				void		readBoundaryInfo(string filename, TPdModel& pdModel)
				{
					multimap<int, vector<int> > Nid_SPC;
					map<int, TBoundaryPrescribedMotion> Nid_PRE_MOTION;
					map<int, TLoadNodePoint> Nid_LOAD_POINT;

					ifstream fin(filename, ios::in);
					if (!fin.is_open())
						return;

					string line;
					int data_type = 0;

					while (!fin.eof())
					{
						getline(fin, line);
						if (line.length() < 1)
							continue;
						if (line.substr(0, 1) == string("$"))
							continue;

						if (string("*BOUNDARY_SPC_NODE") == line)
						{
							streampos pos = fin.tellg();
							while (getline(fin, line))
							{
								if (line.substr(0, 1) == string("*") ||
									line.substr(0, 1) == string("$"))
								{
									fin.seekg(pos);
									break;
								}
								int nid = ParseString<int>(line.substr(0, 10));
								int nid_local = map_fem_nid_global_to_fem_nid_local[nid];

								for (int i = 1; i <= 6; ++i)
								{
									int dof = ParseString<int>(line.substr(10 + 10 * i, 10));
									if (dof == 1)
									{
										pdModel.PdMeshCore().Node(nid_local).BoundarySpcNode().push_back(TBoundarySpcNode(i));
									}
								}
								
								pos = fin.tellg();
							}

							continue;
						}
						else if (string("*BOUNDARY_PRESCRIBED_MOTION_NODE") == line)
						{
							streampos pos = fin.tellg();

							while (getline(fin, line))
							{
								if (line.substr(0, 1) == string("*") ||
									line.substr(0, 1) == string("$"))
								{
									fin.seekg(pos);
									break;
								}
								int nid = ParseString<int>(line.substr(0, 10));
								int nid_local = map_fem_nid_global_to_fem_nid_local[nid];

								TBoundaryPrescribedMotion bpm;
								bpm.Dof() = ParseString<int>(line.substr(10, 10));
								bpm.Vda() = ParseString<int>(line.substr(20, 10));
								bpm.Lcid() = ParseString<int>(line.substr(30, 10));
								bpm.Sf() = ParseString<double>(line.substr(40, 10));
								bpm.MotionType() = "NODE";

								pdModel.PdMeshCore().Node(nid_local).BoundaryPreMotion().push_back(bpm);

								pos = fin.tellg();
							}

							continue;
						}
						else if (string("*LOAD_NODE_POINT") == line)
						{
							streampos pos = fin.tellg();
							while (getline(fin, line))
							{
								if (line.substr(0, 1) == string("*") ||
									line.substr(0, 1) == string("$"))
								{
									fin.seekg(pos);
									break;
								}
								int nid = ParseString<int>(line.substr(0, 10));
								int nid_local = map_fem_nid_global_to_fem_nid_local[nid];
								int dof = ParseString<int>(line.substr(10, 10));
								/************************************************************************/
								/* 5 6 7 分别代表X Y Z的弯矩                                  */
								/************************************************************************/
								if (dof > 4)
								{
									dof -= 1;
								}
								int lcid = ParseString<int>(line.substr(20, 10));
								double sf = ParseString<double>(line.substr(30, 10));

								TLoadNodePoint loadnodepoint;
								loadnodepoint.Dof() = dof;
								loadnodepoint.Lcid() = lcid;
								loadnodepoint.Sf() = sf;

								pdModel.PdMeshCore().Node(nid_local).LoadNodePoint().push_back(loadnodepoint);

								pos = fin.tellg();
							}

							continue;
						}						
					}

					fin.close();
				}
				void		readCalculateInfo(string filename, TPdModel& pdModel)
				{
					ifstream fin(filename, ios::in);
					if (!fin.is_open())
						return;

					string line;
					int data_type = 0;

					while (!fin.eof())
					{
						getline(fin, line);
						if (line.length() < 1)
							continue;
						if (line.substr(0, 1) == string("$"))
							continue;
						
						if (string("*MAT_ELASTIC") == line)
						{
							getline(fin, line);
							int mid = ParseString<int>(line.substr(0, 10));
							double density = ParseString<double>(line.substr(10, 10));
							double elastic_modulus = ParseString<double>(line.substr(20, 10));
							double poisson_ratio = ParseString<double>(line.substr(30, 10));
							double stress_tensile = ParseString<double>(line.substr(40, 10));
							double stress_compressive = ParseString<double>(line.substr(50, 10));
							TMaterial mat;
							mat.Id() = mid;
							mat.Name() = "MAT_ELASTIC";
							mat.InsertMatData("Rho", density);
							mat.InsertMatData("E", elastic_modulus); 
							mat.InsertMatData("PR", poisson_ratio);
							mat.InsertMatData("STRESS_TENSILE", stress_tensile);
							mat.InsertMatData("STRESS_COMPRESSIVE", stress_compressive);
							pdModel.AddMaterial(mat);

							continue;
						}
						else if (string("*MAT_ORTHOTROPIC_ELASTIC") == line)
						{
							getline(fin, line);
							int mid = ParseString<int>(line.substr(0, 10));
							double Rho = ParseString<double>(line.substr(10, 10));
							double EA = ParseString<double>(line.substr(20, 10));
							double EB = ParseString<double>(line.substr(30, 10));
							double EC = ParseString<double>(line.substr(40, 10));
							double PRBA = ParseString<double>(line.substr(50, 10));
							double PRBB = ParseString<double>(line.substr(60, 10));
							double PRBC = ParseString<double>(line.substr(70, 10));
							getline(fin, line);
							double GAB = ParseString<double>(line.substr(0, 10));
							double GBC = ParseString<double>(line.substr(0, 10));
							double GCA = ParseString<double>(line.substr(0, 10));
							getline(fin, line);
							getline(fin, line);
							double BETA = ParseString<double>(line.substr(60, 10));

							TMaterial mat;
							mat.Id() = mid;
							mat.Name() = "MAT_ORTHOTROPIC_ELASTIC";
							mat.InsertMatData("Rho", Rho);
							mat.InsertMatData("EA", EA);
							mat.InsertMatData("EB", EB);
							mat.InsertMatData("EC", EC);
							mat.InsertMatData("PRBA", PRBA);
							mat.InsertMatData("PRBB", PRBB);
							mat.InsertMatData("PRBC", PRBC);
							mat.InsertMatData("GAB", GAB);
							mat.InsertMatData("GBC", GBC);
							mat.InsertMatData("GCA", GCA);
							mat.InsertMatData("BETA", BETA);
							pdModel.AddMaterial(mat);

							continue;
						}
						else if (string("*MAT_RIGID") == line)
						{
							getline(fin, line);
							int mid = ParseString<int>(line.substr(0, 10));
							double Rho = ParseString<double>(line.substr(10, 10));
							double E = ParseString<double>(line.substr(20, 10));
							double PR = ParseString<double>(line.substr(30, 10));
							getline(fin, line);
							double CMO = ParseString<double>(line.substr(0, 10));
							double CON1 = ParseString<double>(line.substr(10, 10));
							double CON2 = ParseString<double>(line.substr(20, 10));
							TMaterial mat;
							mat.Id() = mid;
							mat.Name() = "MAT_RIGID";
							mat.InsertMatData("Rho", Rho);
							mat.InsertMatData("E", E);
							mat.InsertMatData("PR", PR);
							mat.InsertMatData("CMO", CMO);
							mat.InsertMatData("CON1", CON1);
							mat.InsertMatData("CON2", CON2);
							pdModel.AddMaterial(mat);

							continue;
						}
						else if (string("*PART") == line)
						{
							streampos pos = fin.tellg();
							while (getline(fin, line))
							{
								if (line.substr(0, 1) == string("*") ||
									line.substr(0, 1) == string("$"))
								{
									fin.seekg(pos);
									break;
								}
							string partname = line;		//	PART_NAME
							getline(fin, line);
							int partid = ParseString<int>(line.substr(0, 10));
							int propid = ParseString<int>(line.substr(10, 10));
							int matid = ParseString<int>(line.substr(20, 10));
							TPart part;
							part.Id() = partid;
							part.MaterialId() = matid;
							part.SectionId() = propid;
							pdModel.AddPart(part);
							}

							continue;
						}
						else if (string("*SECTION_SHELL") == line)
						{
							getline(fin, line);
							int pid = ParseString<int>(line.substr(0, 10));		//	SECTION_ID
							int ELEFORM = ParseString<int>(line.substr(10, 10));
							getline(fin, line);
							double THICKNESS = ParseString<double>(line.substr(0, 10));
							TSection sec(pid);
							if (ELEFORM == 13)
							{
								sec.Type() = string("PLANE_STRAIN");
							}
							else
							{
								sec.Type() = string("PLANE_STRESS");
							}
							sec.InsertSectionData("THICKNESS", THICKNESS);
							pdModel.AddSection(sec);

							continue;
						}
						else if (string("*DEFINE_CURVE") == line)
						{
							getline(fin, line);
							int cid = ParseString<int>(line.substr(0, 10));
							TCurve curve(cid);

							streampos pos = fin.tellg();
							while (getline(fin, line))
							{
								if (line.substr(0, 1) == string("*"))
								{
									fin.seekg(pos);
									break;
								}
								double x = ParseString<double>(line.substr(0, 20));
								double y = ParseString<double>(line.substr(20, 20));
								curve.AddPoint(x, y);

								pos = fin.tellg();
							}

							pdModel.AddCurve(curve);

							continue;
						}
						else if (string("*ELEMENT_SEATBELT") == line)
						{
							streampos pos = fin.tellg();
							while (getline(fin, line))
							{
								if (line.substr(0, 1) == string("*"))
								{
									fin.seekg(pos);
									break;
								}
								int n1 = ParseString<int>(line.substr(16, 8));
								int n2 = ParseString<int>(line.substr(24, 8));

								Vector3d p1 = pdModel.PdMeshCore().Node(map_fem_nid_global_to_fem_nid_local[n1]).Coordinate();
								Vector3d p2 = pdModel.PdMeshCore().Node(map_fem_nid_global_to_fem_nid_local[n2]).Coordinate();

								TCrevice crevice(p1, p2);
								pdModel.AddCrevice(crevice);

								pos = fin.tellg();
							}

							continue;
						}										
						else if (string("*SET_SHELL_LIST_TITLE") == line)
						{
							getline(fin, line);
							string setName;
							stringstream input(line);
							input >> setName;

							getline(fin, line);
							int setId = ParseString<int>(line.substr(0, line.length()-1));

							streampos pos = fin.tellg();
							while (getline(fin, line))
							{
								if (line.substr(0, 1) == string("*") ||
									line.substr(0, 1) == string("$"))
								{
									fin.seekg(pos);
									break;
								}

								int num = (int)(line.length()) / 10;
								for (int i = 0; i < num; ++i)
								{
									int eid = ParseString<int>(line.substr(0 + 10 * i, 10));
									if (eid != 0)
									{
										int eid_local = map_fem_eid_global_to_fem_eid_local[eid];

										if (setName == "PD")
										{
											pdModel.PdMeshCore().Element(eid_local).Alpha() = 1.0;
										} 
									}
								}

								pos = fin.tellg();
							}

							continue;
						}
					}
					fin.close();
				}
			private:
				template<typename T>
				T			ParseString(string str)
				{
					stringstream sstr;
					sstr << str;
					T Ret;
					if (str == "          " ||
						str == "        ")
					{
						Ret = 0;
					}
					else
					{
						sstr >> Ret;
					}
					return Ret;
				}
			private:
				map<int, int>		map_fem_nid_global_to_fem_nid_local;		//	FEM NODE GLOBAL    <-> FEM NODE LOCAL
				map<int, int>		map_fem_eid_global_to_fem_eid_local;		//	FEM ELEMENT GLOBAL <-> FEM ELEMENT LOCAL
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif