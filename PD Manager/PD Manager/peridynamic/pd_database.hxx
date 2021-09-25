/*All rights reserved

== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
--Please append file description informations here --

== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
Date            Name                    Description of Change

$HISTORY$
== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
*/

#ifndef DLUT_SAE_PERIDYNAMIC_PD_DATABASE_HXX_20171220
#define DLUT_SAE_PERIDYNAMIC_PD_DATABASE_HXX_20171220

#include <map>
#include <string>

#include "UMF/umfSolver.h"

#include "pd_base_toolkit.hxx"
#include "fem_shape_functions.hxx"
#include "fem_database.hxx"
#include <ppl.h>
using namespace concurrency;

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			static double RATIO_OF_HORIZON_MESHSIZE = 3.0;
			static double HORIZON = 3.0;
			static bool USE_CONSTANT_HORIZON = TRUE;
			static int MAX_INTERATOR_NUMS = 100;
			static double CONVERGENCE_FACTOR = 0.001;
			/************************************************************************/
			/* 采用Beam PD时，节点采用6个自由度,采用BBPD时，节点采用3个自由度       */
			/************************************************************************/
			class TPdCalculateParas
			{
			public:
				TPdCalculateParas()
				{
					D_PD.setZero();

					sed_criterion = 0;
				}
				~TPdCalculateParas() {}
			public:
				Eigen::Matrix4d D_PD;
				
				double	sed_criterion;
			};

			typedef TIntegrationPointTemp<TStrain, TStress> TBondPoint;
			typedef TIntegrationPointTemp<TStrain_Bond, TStress_Bond> TBondIntegrationPoint;
			class TPdBond
			{
			public:
				TPdBond(TBondPoint& xi, TBondPoint& xj, const Matrix3d& bond_local_coor_sys)
					: m_xi(xi), m_xj(xj)
				{
					m_local_coord_system = bond_local_coor_sys;
					m_IP.clear();
					m_IP.resize(IP_COUNT_1D);
					m_micro_potential = 0;
					m_b_is_valid = true;
				}
			public:
				TBondPoint&			Xi() { return m_xi; }
				const TBondPoint&	Xi() const { return m_xi; }
				TBondPoint&			Xj() { return m_xj; }
				const TBondPoint&	Xj() const { return m_xj; }
				double				BondLength() const
				{
					return Distance_2pt(m_xi.Coordinate(), m_xj.Coordinate());
				}
			public:
				Matrix3d&			LocalCoorSystem() { return m_local_coord_system; }
				const Matrix3d&		LocalCoorSystem() const { return m_local_coord_system; }
				//	键的局部坐标系与全局坐标系的转换矩阵
				Eigen::MatrixXd		T_b() const
				{
					Eigen::MatrixXd Res;
					Res.resize(12, 12);
					Res.setZero();

					Res.block(0, 0, 3, 3) = m_local_coord_system;
					Res.block(3, 3, 3, 3) = m_local_coord_system;
					Res.block(6, 6, 3, 3) = m_local_coord_system;
					Res.block(9, 9, 3, 3) = m_local_coord_system;

					return Res;
				}
			public:
				TBondIntegrationPoint&			IP(int index) { return m_IP[index]; }
				const TBondIntegrationPoint&	IP(int index) const { return m_IP[index]; }
			public:
				double						MicroPotential() const { return m_micro_potential; }
				double&						MicroPotential() { return m_micro_potential; }
				bool						IsValid() const { return m_b_is_valid; }
				void						MakeFailure() { m_b_is_valid = false; }
			private:
				double						m_micro_potential;
				bool						m_b_is_valid;
			private:
				Matrix3d					m_local_coord_system;	//	Local coordinate system of the bond
			private:
				vector<TBondIntegrationPoint>	m_IP;					//	Integration Points of this bond
			private:
				//	组成TPdBond的是两个单元中的积分点
				TBondPoint& m_xi;					//	Xi;
				TBondPoint& m_xj;					//	Xj;
			};

			class TPdFamilyElement
			{
			public:
				TPdFamilyElement(int eid = -1, double volume_index = 0)
				{
					m_eid_j = eid;
					m_volume_index = volume_index;
					m_bonds.clear();
					m_b_update = true;
				}
			public:
				int						Id() const { return m_eid_j; }
				double&					VolumeIndex() { return m_volume_index; }
				double					VolumeIndex() const { return m_volume_index; }
			public:
				void					AddBond(TBondPoint& xi, TBondPoint& xj, const Matrix3d& bond_local_coor_sys = Matrix3d().Identity())
				{
					//	两个积分点不为同一个值才能形成一根bond
					if (&xi != &xj)
					{
						m_bonds.push_back(TPdBond(xi, xj, bond_local_coor_sys));
					}
				}
				TPdBond&				Bond(int index) { return m_bonds[index]; }
				const TPdBond&			Bond(int index) const { return m_bonds[index]; }
				vector<TPdBond>&		Bonds() { return m_bonds; }
				const vector<TPdBond>&	Bonds() const { return m_bonds; }
			public:
				SingleStiffness&		SK() { return m_single_stiffness; }
				const SingleStiffness&	SK() const { return m_single_stiffness; }
				Eigen::VectorXd&		ForceOfBond() { return m_force; }
				const Eigen::VectorXd&	ForceOfBond() const { return m_force; }
				const bool				ShouldBeUpdate() const { return m_b_update; }
				bool&					ShouldBeUpdate() { return m_b_update; }
			private:
				int						m_eid_j;				//	Element J
				double					m_volume_index;				//	Modified volume
				vector<TPdBond>			m_bonds;				//	TPdBonds.
				SingleStiffness			m_single_stiffness;		//	Single stiffness matrix of Bond_ij
				Eigen::VectorXd			m_force;				//	FORCE vector of Element_ij
				bool					m_b_update;				//	SK should be updated?
			};

			//typedef map<int, TPdFamilyElement>	MAP_NJ_FAMILY_ELEMENT;
			//typedef pair<int, TPdFamilyElement> PAIR_NJ_FAMILY_ELEMENT;

			typedef vector<TPdFamilyElement> VEC_FAMILY_ELEMENT;

			typedef TNodeBase TPdNode;
			//	Element
			class TPdElement : public TElementBase
			{
			public:
				TPdElement(vector<TNodeBase>& vecNode) : TElementBase(vecNode)
				{
					m_vec_family_elements.clear();
					m_d_init_bond_nums = 0;
					m_alpha = 0;
				}
			public:
				void					Dispose()
				{
					TElementBase::Dispose();
				}
			/*	const TPdElement& operator=(const TPdElement& right)
				{
					TElementBase::operator=(right);
					m_map_family_elements = right.m_map_family_elements;
					m_pd_paras = right.m_pd_paras;
					m_d_init_bond_nums = right.m_d_init_bond_nums;

					return *this;
				}*/
			public:
				//	The family nodes informations and operations
	/*			int						FamilyElementCount() const { return (int)m_map_family_elements.size(); }
				MAP_NJ_FAMILY_ELEMENT&	FamilyElements() { return m_map_family_elements; }
				const MAP_NJ_FAMILY_ELEMENT& FamilyElements() const { return m_map_family_elements; }
				TPdFamilyElement&		FamilyElement(int eId) { return m_map_family_elements[eId]; }
				void					InsertFamilyElement(int eId, double volume_index) { m_map_family_elements.insert(pair<int, TPdFamilyElement>(eId, TPdFamilyElement(eId, volume_index))); }
				void					DeleteFamilyElement(int eId) { m_map_family_elements.erase(eId); }
				void					ClearFamilyElements() { m_map_family_elements.clear(); }*/

				int						FamilyElementCount() const { return (int)m_vec_family_elements.size(); }
				VEC_FAMILY_ELEMENT&		FamilyElements() { return m_vec_family_elements; }
				const VEC_FAMILY_ELEMENT& FamilyElements() const { return m_vec_family_elements; }
				void					InsertFamilyElement(int eId, double volume_index) { m_vec_family_elements.push_back(TPdFamilyElement(eId, volume_index)); }
				void					ClearFamilyElements() { m_vec_family_elements.clear(); }
				TPdFamilyElement&		LastOfFamilyElement() { return m_vec_family_elements.back(); }
			public:
				const TPdCalculateParas&	CalParas() const { return m_pd_paras; }
				TPdCalculateParas&			CalParas() { return m_pd_paras; }
			public:
				void					InitDamageIndex()
				{
					m_d_init_bond_nums = 0;
					const VEC_FAMILY_ELEMENT& familyElems = FamilyElements();
					for (const TPdFamilyElement& family_elem : familyElems)
					{
						const vector<TPdBond>& bonds = family_elem.Bonds();
						for (const TPdBond& bond : bonds)
						{
							if (bond.IsValid())
							{
								m_d_init_bond_nums += 1;
							}
						}
					}
				}
				double					DamageIndex() const
				{
					double invalid_bond_nums = 0;
					const VEC_FAMILY_ELEMENT& familyElems = FamilyElements();
					for (const TPdFamilyElement& family_elem : familyElems)
					{
						const vector<TPdBond>& bonds = family_elem.Bonds();
						for (const TPdBond& bond : bonds)
						{
							if (!bond.IsValid())
							{
								invalid_bond_nums += 1;
							}
						}
					}
					return invalid_bond_nums / m_d_init_bond_nums;
				}

			public:
				double					Alpha() const { return m_alpha; }
				double&					Alpha() { return m_alpha; }

			private:
//				MAP_NJ_FAMILY_ELEMENT	m_map_family_elements;				//	Bond informations
				VEC_FAMILY_ELEMENT		m_vec_family_elements;  			//	Bond informations
				TPdCalculateParas		m_pd_paras;							//	Calculate parameters of PD Element
				double					m_d_init_bond_nums;					//	Initial volumes for damage calculation

			private:
				double					m_alpha;							//	Alpha index for the coupling model
			};
			
			typedef TMeshCoreTemplate<TPdNode, TPdElement> TPdMeshCore;

			class TPdDataCollector
			{				
			public:
				TPdDataCollector() 
				{ 
					Initialize();
				}
				~TPdDataCollector() {}
			public:
				void						Initialize()
				{
					m_pd_meshcore.Initialize();

					m_vec_parts.clear();
					m_vec_materials.clear();
					m_vec_sections.clear();
					m_vec_curves.clear();

					m_vec_crevice.clear();
				}
			public:
				TPdMeshCore&				PdMeshCore() { return m_pd_meshcore; }
				const TPdMeshCore&			PdMeshCore() const { return m_pd_meshcore; }
			public:
				TPart&						Part(int part_id)
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_parts.size()); ++loop)
					{
						if (part_id == m_vec_parts[loop].Id())
						{
							cur_loop = loop;
							break;
						}
					}
					return m_vec_parts[cur_loop];
				}
				const TPart&				Part(int part_id) const
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_parts.size()); ++loop)
					{
						if (part_id == m_vec_parts[loop].Id())
						{
							cur_loop = loop;
							break;
						}
					}
					return m_vec_parts[cur_loop];
				}
				void						AddPart(const TPart& part)
				{
					bool b_exist = false;
					for (int loop = 0; loop < (int)(m_vec_parts.size()); ++loop)
					{
						TPart& cur_part = m_vec_parts[loop];
						//	如果有同样ID号的PART，则更新信息即可
						if (cur_part.Id() == part.Id())
						{
							m_vec_parts[loop] = part;
							b_exist = true;
							break;
						}
					}
					//	如果不存在同样ID的PART，则新插入一个
					if (!b_exist)
					{
						m_vec_parts.push_back(part);
					}
				}
				int							PartCounts() const { return (int)(m_vec_parts.size()); }
				int							PartId(int i_count) const { return m_vec_parts[i_count].Id(); }
				const vector<TPart>&		Parts() { return m_vec_parts; }

				TMaterial&					Material(int mat_id)
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_materials.size()); ++loop)
					{
						if (mat_id == m_vec_materials[loop].Id())
						{
							cur_loop = loop;
						}
					}
					return m_vec_materials[cur_loop];
				}
				const TMaterial&			Material(int mat_id) const
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_materials.size()); ++loop)
					{
						if (mat_id == m_vec_materials[loop].Id())
						{
							cur_loop = loop;
						}
					}
					return m_vec_materials[cur_loop];
				}
				bool						MaterialExist(int mat_id)
				{
					bool res = false;
					for (int loop = 0; loop < (int)(m_vec_materials.size()); ++loop)
					{
						if (mat_id == m_vec_materials[loop].Id())
						{
							res = true;
							break;
						}
					}
					return res;
				}
				void						AddMaterial(const TMaterial& mat)
				{
					m_vec_materials.push_back(mat);
				}
				int							MaterialCounts() const { return (int)(m_vec_materials.size()); }
				int							MaterialId(int i_count) const { return m_vec_materials[i_count].Id(); }
				const vector<TMaterial>&	Materials() { return m_vec_materials; }

				TSection&					Section(int sec_id)
				{
					int cur_loop = -1;
					for (int loop = 0; loop < (int)(m_vec_sections.size()); ++loop)
					{
						if (sec_id == m_vec_sections[loop].Id())
						{
							cur_loop = loop;
						}
					}
					assert(cur_loop >= 0);
					return m_vec_sections[cur_loop];
				}
				const TSection&				Section(int sec_id) const
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_sections.size()); ++loop)
					{
						if (sec_id == m_vec_sections[loop].Id())
						{
							cur_loop = loop;
						}
					}
					return m_vec_sections[cur_loop];
				}
				bool						SectionExist(int sec_id)
				{
					bool res = false;
					for (int loop = 0; loop < (int)(m_vec_sections.size()); ++loop)
					{
						if (sec_id == m_vec_sections[loop].Id())
						{
							res = true;
							break;
						}
					}
					return res;
				}
				void						AddSection(const TSection& sec)
				{
					m_vec_sections.push_back(sec);
				}
				int							SectionCounts() const { return (int)(m_vec_sections.size()); }
				int							SectionId(int i_count) const { return m_vec_sections[i_count].Id(); }
				const vector<TSection>&		Sections() { return m_vec_sections; }
				
				TCurve&						Curve(int cur_id)
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)m_vec_curves.size(); ++loop)
					{
						TCurve& curve = m_vec_curves[loop];
						if (cur_id == curve.Id())
						{
							cur_loop = loop;
							break;
						}
					}
					return m_vec_curves[cur_loop];
				}
				const TCurve&				Curve(int cur_id) const
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)m_vec_curves.size(); ++loop)
					{
						const TCurve& curve = m_vec_curves[loop];
						if (cur_id == curve.Id())
						{
							cur_loop = loop;
							break;
						}
					}
					return m_vec_curves[cur_loop];
				}
				bool						CurveExist(int cur_id)
				{
					bool res = false;
					for (int loop = 0; loop < (int)m_vec_curves.size(); ++loop)
					{
						const TCurve& curve = m_vec_curves[loop];
						if (cur_id == curve.Id())
						{
							res = true;
							break;
						}
					}
					return res;
				}
				void						AddCurve(const TCurve& cur)
				{
					m_vec_curves.push_back(cur);
				}
				int							CurveCounts() const { return (int)(m_vec_curves.size()); }
				int							CurveId(int i_count) { return m_vec_curves[i_count].Id(); }
				const vector<TCurve>&		Curves() { return m_vec_curves; }
				
				//	Add the initial crevice into the PD model, and the crevice is *ELEMENT_SEATBELT
				void						AddCrevice(const TCrevice& cre)
				{
					m_vec_crevice.push_back(cre);
				}
				const vector<TCrevice>&		Crevice() const { return m_vec_crevice; }
			protected:
				TPdMeshCore							m_pd_meshcore;			//	PD MESH CORE

				vector<TPart>						m_vec_parts;			//	*PART
				vector<TMaterial>					m_vec_materials;		//	*MAT
				vector<TSection>					m_vec_sections;			//	*SECTION
				vector<TCurve>						m_vec_curves;			//	*CURVE
				vector<TCrevice>					m_vec_crevice;			//	*ELEMENT_SEATBELT
			};

			//	PD instantiation model
			class TPdModel : public TPdDataCollector
			{
			public:
				//	Update the *PART informations into PD MODEL
				void				UpdatePartInfo()
				{	
					double start, end, cost;
					start = clock();
					if (PartCounts() == 0)
					{
						cout << "ERROR: Have no any part information in this LSDYNA file!" << endl;
						return;
					}
					
					//	获取Part的信息：单元编号、Section信息、材料信息等
					for (int pid = 0; pid < PartCounts(); ++pid)
					{
						TPart& part = Part(PartId(pid));
						
						set<int> eids = m_pd_meshcore.GetElementIdsByPart(part.Id());
						part.AddElementId(eids);
						//	设置单元的AREA
						const TSection& section = Section(part.SectionId());
						string stype = section.Type();
						if (stype == "PLANE_STRESS" ||
							stype == "PLANE_STRAIN")
						{
							double thickness = section.GetSectionValue("THICKNESS");
							for (int eid : eids)
							{
								TPdElement& pd_element = m_pd_meshcore.Element(eid);
								pd_element.Thickness() = thickness;
							}
						}
					}
					
					//	更新材料信息
					for (const TPart& part : Parts())
					{
						if (!MaterialExist(part.MaterialId()))
						{
							cout << "Material " << part.MaterialId() << " is not exist in this KEYWORD file." << endl;
							return;
						}
						const TMaterial& material = Material(part.MaterialId());
						if (material.Name() == string("MAT_RIGID"))
						{
							int cmo = (int)(material.GetMatValue("CMO"));
							if (cmo == 1)
							{
								int con1 = (int)(material.GetMatValue("CON1"));
								int con2 = (int)(material.GetMatValue("CON2"));
								set<int> nids = part.GetElementIds();
								for (int nid : nids)
								{
									PdMeshCore().Node(nid).BoundarySpcNode().push_back(TBoundarySpcNode(1));
									/*
									switch (con1)
									{
									case 1:
										AddBounarySpcNode(TBoundarySpcNode(nid, 1));
										break;
									case 2:
										AddBounarySpcNode(TBoundarySpcNode(nid, 2));
										break;
									case 3:
										AddBounarySpcNode(TBoundarySpcNode(nid, 3));
										break;
									case 4:
										AddBounarySpcNode(TBoundarySpcNode(nid, 1));
										AddBounarySpcNode(TBoundarySpcNode(nid, 2));
										break;
									case 5:
										AddBounarySpcNode(TBoundarySpcNode(nid, 2));
										AddBounarySpcNode(TBoundarySpcNode(nid, 3));
										break;
									case 6:
										AddBounarySpcNode(TBoundarySpcNode(nid, 3));
										AddBounarySpcNode(TBoundarySpcNode(nid, 1));
										break;
									case 7:
										AddBounarySpcNode(TBoundarySpcNode(nid, 1));
										AddBounarySpcNode(TBoundarySpcNode(nid, 2));
										AddBounarySpcNode(TBoundarySpcNode(nid, 3));
										break;
									default:
										break;
									}*/
								}
							}
						}
					}
					end = clock();
					cost = end - start;
					cout << "Update Part Informations, Cost time = " << cost / 1000 << endl;
				}
				void				UpdateFamilyInParts()
				{
					double start, end, cost;
					start = clock();
					//	First, delete all exist family information
					for (int loop = 0; loop < m_pd_meshcore.ElementCount(); ++loop)
					{
						m_pd_meshcore.Element(loop).ClearFamilyElements();
					}

					//	Then, rebuild the family information for every node
					for (const TPart& part : Parts())
					{
						int sid = part.SectionId();
						const TSection& section = Section(sid);
						string stype = section.Type();

						const set<int>& eleIds = part.GetElementIds();
						for (int ei : eleIds)
						{
							TPdElement& element_i = m_pd_meshcore.Element(ei);
							TCoordinate coor_i = element_i.CoordinateInElement(0, 0);

							double dx = element_i.SideLength();
							//	Extend the search range to include all nodes

							double dis_for_juge = 0;
							double horizon = 0;
							if (DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON)
							{
								dis_for_juge = DLUT::SAE::PERIDYNAMIC::HORIZON + dx + ERR_VALUE;
								horizon = DLUT::SAE::PERIDYNAMIC::HORIZON;
							} 
							else
							{
								dis_for_juge = (DLUT::SAE::PERIDYNAMIC::RATIO_OF_HORIZON_MESHSIZE + 1.0) * dx + ERR_VALUE;
								horizon = DLUT::SAE::PERIDYNAMIC::RATIO_OF_HORIZON_MESHSIZE * dx;
							}

							int iter_for_adj = (int)(ceil(horizon / dx));
							set<int> ejIds = m_pd_meshcore.GetAdjElements(ei, iter_for_adj);
							for (int ej : ejIds)
							{								
								TPdElement& element_j = m_pd_meshcore.Element(ej);
								TCoordinate coor_j = element_j.CoordinateInElement(0, 0);

								if ((abs(coor_i.x() - coor_j.x()) < dis_for_juge) &&
									(abs(coor_i.y() - coor_j.y()) < dis_for_juge) &&
									(abs(coor_i.z() - coor_j.z()) < dis_for_juge))
								{
									double Len = Distance_2pt(coor_i, coor_j);
									if (Len < dis_for_juge)
									{
										Vector3d cij;
										double intersect_volume_index = CalModifiedVolume(element_i, element_j, cij);
										if (intersect_volume_index > 0)
										{
											element_i.InsertFamilyElement(ej, intersect_volume_index);
											for (int is = 0; is < IP_COUNT_2D; ++is)
											{
												TIntegrationPoint& xi = element_i.IP(is);
												for (int js = 0; js < IP_COUNT_2D; ++js)
												{
													TIntegrationPoint& xj = element_j.IP(js);
													if (element_i.Id() == element_j.Id() && xi.Index() == xj.Index())
													{
														continue;
													}
													Matrix3d bond_local_coor_sys;
													bond_local_coor_sys.setZero();
													Vector3d lx = xj.Coordinate() - xi.Coordinate();
													Vector3d lz = element_i.LocalCoorSystem().block(2, 0, 1, 3).transpose();
													Vector3d ly = Fork_Multi<Vector3d, Vector3d>(lz, lx);
													Normalizer<Vector3d>(lx);
													Normalizer<Vector3d>(ly);
													Normalizer<Vector3d>(lz);
													bond_local_coor_sys.block(0, 0, 1, 3) = lx.transpose();
													bond_local_coor_sys.block(1, 0, 1, 3) = ly.transpose();
													bond_local_coor_sys.block(2, 0, 1, 3) = lz.transpose();

													//	刚插入的FamilyElement进行Bond信息更新
													element_i.LastOfFamilyElement().AddBond(xi, xj, bond_local_coor_sys);
												}
											}											
											//	根据单元位置信息，计算耦合标量函数的信息alpha										
											if (element_i.Alpha() > (1.0 - ERR_VALUE))
											{
												element_i.Alpha() = 1;
												double alpha = CalAlpha(Len, horizon);
												if (element_j.Alpha() < alpha)
												{
													element_j.Alpha() = alpha;
												}
											}
										}
									}
								}
							}
						}
						//	Thirdly, divide the total domain into PD, FEM, and MORPHING domain, respectively.
						//for (int ei : eleIds)
						//{
						//	TPdElement& element_i = m_pd_meshcore.Element(ei);
						//	MAP_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
						//	if (element_i.Alpha() > (1.0 - ERR_VALUE))
						//	{
						//		bool is_pd = true;
						//		for (MAP_NJ_FAMILY_ELEMENT::iterator iter = familyElements.begin();
						//			iter != familyElements.end(); ++iter)
						//		{
						//			TPdElement& element_j = PdMeshCore().Element((*iter).first);
						//			if (element_j.Alpha() < 1.0)
						//			{
						//				is_pd = false;
						//				break;
						//			}
						//		}
						//		if (is_pd)
						//		{
						//			element_i.AnalysisElementType() = PD_ELEMENT;
						//			//	将所有PD单元设置成离散单元
						//			m_pd_meshcore.AddSeparateElement(ei);
						//		}
						//		else
						//		{
						//			element_i.AnalysisElementType() = MORPHING_ELEMENT;
						//		}
						//	}
						//	else if (element_i.Alpha() < ERR_VALUE)
						//	{
						//		element_i.AnalysisElementType() = FEM_ELEMENT;
						//		bool is_fem = true;
						//		for (MAP_NJ_FAMILY_ELEMENT::iterator iter = familyElements.begin();
						//			iter != familyElements.end(); ++iter)
						//		{
						//			TPdElement& element_j = PdMeshCore().Element((*iter).first);
						//			if (element_j.Alpha() > ERR_VALUE)
						//			{
						//				is_fem = false;
						//				break;
						//			}
						//		}
						//		if (is_fem)
						//		{
						//			element_i.AnalysisElementType() = FEM_ELEMENT;
						//		}
						//		else
						//		{
						//			element_i.AnalysisElementType() = MORPHING_ELEMENT;
						//		}
						//	}
						//	else
						//	{
						//		element_i.AnalysisElementType() = MORPHING_ELEMENT;
						//	}
						//}
					}

					end = clock();
					cost = end - start;
					cout << "Update Family Informations, Cost time = " << cost / 1000 << endl;
				}
			public:
				//	Refresh the *BOUNDARY_SPC_NODE informations into PD MODEL
				void				RefreshBoundarySpcInfo()
				{
					const set<int>& nids = PdMeshCore().GetNodeIdsByAll();
					for (int nid : nids)
					{
						TPdNode& node = PdMeshCore().Node(nid);
						const vector<TBoundarySpcNode>& BSNs = node.BoundarySpcNode();
						for (const TBoundarySpcNode& tbsn : BSNs)
						{
							int dof = tbsn.Dof() - 1;
							node.Displacement()(dof) = 0;
						}
					}
				}
				//	Refresh the *INITIAL_VELOCITY informations into PD MODEL
				void				RefreshInitVelocityInfo()
				{
					const set<int>& nids = PdMeshCore().GetNodeIdsByAll();
					for (int nid : nids)
					{
						TPdNode& node = PdMeshCore().Node(nid);
						const vector<TInitialVelocityNode>& IVNs = node.InitVelocity();
						for (const TInitialVelocityNode& tivn : IVNs)
						{
							node.Velocity() = tivn.Velocity();
						}
					}					
				}
				//	Refresh the *LOAD_NODE_POINT informations into PD MODEL
				void				RefreshLoadNodeInfo(double cur_time)
				{				
					const set<int>& nids = PdMeshCore().GetNodeIdsByAll();
					for (int nid : nids)
					{
						TPdNode& node = PdMeshCore().Node(nid);
						const vector<TLoadNodePoint>& LNPs = node.LoadNodePoint();
						for (const TLoadNodePoint& tlnp : LNPs)
						{
							int curid = tlnp.Lcid();

							double value_force = 0;
							if (CurveExist(curid))
							{
								TCurve& curve = Curve(curid);
								value_force = curve.GetValueByX(cur_time) * tlnp.Sf();
							}
							else
							{
								value_force = tlnp.Sf();
							}

							switch (tlnp.Dof())
							{
							case 1:
							{
								node.OuterForce().x() = value_force;
								break;
							}
							case 2:
							{
								node.OuterForce().y() = value_force;
								break;
							}
							case 3:
							{
								node.OuterForce().z() = value_force;
								break;
							}
							default:
								break;
							}
						}
					}
				}
				//	Refresh the *BOUNDARY_PRESCRIBED_MOTION_NODE information into the PD MODEL 
				void				RefreshPreMotionInfo(double cur_time)
				{	
					const set<int>& nids = PdMeshCore().GetNodeIdsByAll();
					for (int nid : nids)
					{
						TPdNode& node = PdMeshCore().Node(nid);
						const vector<TBoundaryPrescribedMotion>& BPMs = node.BoundaryPreMotion();
						for (const TBoundaryPrescribedMotion& bpm : BPMs)
						{
							int curid = bpm.Lcid();
							double value = 0;
							if (CurveExist(curid))
							{
								TCurve& curve = Curve(curid);
								value = curve.GetValueByX(cur_time) * bpm.Sf();
							}
							else
							{
								value = bpm.Sf();
							}

							switch (bpm.Dof())
							{
								//	DOF=1 -> X translation
							case 1:
							{
								switch (bpm.Vda())
								{
									//	VDA = 0 -> VELOCITY
								case 0:
								{
									m_pd_meshcore.Node(nid).Displacement().x() = value * cur_time;
									break;
								}
								//	VDA = 1 -> ACCELERATION
								case 1:
								{
									m_pd_meshcore.Node(nid).Displacement().x() = 0.5 * value * cur_time * cur_time;
								}
								//	VDA = 2 -> DISPLACEMENT
								case 2:
								{
									m_pd_meshcore.Node(nid).Displacement().x() = value;
								}
								default:
									break;
								}
								break;
							}
							//	DOF=2 -> Y translation
							case 2:
							{
								switch (bpm.Vda())
								{
								case 0:
								{
									m_pd_meshcore.Node(nid).Displacement().y() = value * cur_time;
									break;
								}
								case 1:
								{
									m_pd_meshcore.Node(nid).Displacement().y() = 0.5 * value * cur_time * cur_time;
								}
								case 2:
								{
									m_pd_meshcore.Node(nid).Displacement().y() = value;
								}
								default:
									break;
								}
								break;
							}
							//	DOF=3 -> Z translation
							case 3:
							{
								switch (bpm.Vda())
								{
								case 0:
								{
									m_pd_meshcore.Node(nid).Displacement().z() = value * cur_time;
									break;
								}
								case 1:
								{
									m_pd_meshcore.Node(nid).Acceleration().z() = 0.5 * value * cur_time * cur_time;
								}
								case 2:
								{
									m_pd_meshcore.Node(nid).Displacement().z() = value;
								}
								default:
									break;
								}
								break;
							}
							default:
								break;
							}
						}
					}
				}
			public:
				//	Generate initial crevice for PD model
				void				GenerateInitCrevice()
				{
					double start = clock();
				
					for (const TCrevice& crevice : m_vec_crevice)
					{
						for (const TPart& part : Parts())
						{
							const set<int> eids = part.GetElementIds();
							parallel_for_each(eids.begin(), eids.end(), [&](int ei) {
								TPdElement& element_i = m_pd_meshcore.Element(ei);
								const TCoordinate& ei_coord = element_i.CoordinateInElement(0, 0);
								VEC_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
								for (TPdFamilyElement& family_elem : familyElements)
								{
									int ej = family_elem.Id();
									TPdElement& element_j = m_pd_meshcore.Element(ej);
									const TCoordinate& ej_coord = element_j.CoordinateInElement(0, 0);

									bool res = IsTowLineIntersect_xy_plane<Vector3d>(ei_coord, ej_coord, crevice.Start(), crevice.End());
									//	如果单元I与J的中心点连线与初始裂缝相交，则所有的Bond均失效
									if (res)
									{
										vector<TPdBond>& bonds = family_elem.Bonds();
										for (TPdBond& bond : bonds)
										{
											bond.MakeFailure();
										}
									}
								}
							});
						}
					}
					
					double total_time = (clock() - start) / 1000;
					cout << "GenerateInitCrevice():\t\t" << total_time << endl;
				}
				//	Generate initial Damage Value for all PD nodes
				void				GenerateInitDamage()
				{
					double start = clock();
					//	Set FamilyNodeCount for Initial Damage Value
					set<int> eleIds = m_pd_meshcore.GetElementIdsByAll();
					for (int eid : eleIds)
					{
						TPdElement& element = m_pd_meshcore.Element(eid);
						element.InitDamageIndex();
					}

					double total_time = (clock() - start) / 1000;
					cout << "GenerateInitDamage():\t\t" << total_time << endl;
				}
			public:
				double				CalModifiedVolume(const TPdElement& element_i, const TPdElement& element_j, Vector3d& center)
				{
					double volume_of_j = 0;
				
					double Ri = 0;
					if (DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON)
					{
						Ri = DLUT::SAE::PERIDYNAMIC::HORIZON;
					}
					else
					{
						Ri = element_i.SideLength() * DLUT::SAE::PERIDYNAMIC::RATIO_OF_HORIZON_MESHSIZE;
					}

					{
						int partId = element_i.PartId();
						const TPart& part = Part(partId);
						int sectionId = part.SectionId();
						const TSection& section = Section(sectionId);
						double thickness = section.GetSectionValue("THICKNESS");

						vector<int> nids = element_j.NodeIds();
						int nCount = (int)(nids.size());
				
						Vector3d pt1 = m_pd_meshcore.Node(nids[0]).Coordinate();
						Vector3d pt2 = m_pd_meshcore.Node(nids[1]).Coordinate();
						Vector3d pt3 = m_pd_meshcore.Node(nids[2]).Coordinate();

						if (element_j.MeshElementType() == TRIANGLE_ELEMENT)
						{
							volume_of_j = Calculate_Intersect_Area_Circle_Triangle(element_i.CoordinateInElement(0, 0), Ri, pt1, pt2, pt3, center) * thickness;
						}
						else
						{
							Vector3d pt4 = m_pd_meshcore.Node(nids[3]).Coordinate();
							volume_of_j = Calculate_Intersect_Area_Circle_Quadrangle(element_i.CoordinateInElement(0, 0), Ri, pt1, pt2, pt3, pt4, center) * thickness;
						}
					}

					return volume_of_j;
				}				
				double				CalAlpha(double xi, double horizon)
				{
					double alpha = 1 - xi / horizon;
					if (alpha < 0)
					{
						alpha = 0;
					}
					if (alpha > 1)
					{
						alpha = 1;
					}
					return alpha;
				}
			};

			/************************************************************************/
			/* 从 vector< map<int, double> > 转换到 稀疏矩阵                         */
			/************************************************************************/
			void TransVecMap2SparseMatrix(const vector< map<int, double> >& vec_map_matrix, SparseMatrix<double>& sparse_matrix)
			{
				const double ERROR_FOR_SM = ERR_VALUE;
				int nCount = (int)(vec_map_matrix.size());
				sparse_matrix.resize(nCount, nCount);

				vector< Triplet<double> > tri;
				tri.clear();
				for (int i = 0; i < (int)(vec_map_matrix.size()); ++i)
				{
					for (const pair<int, double>& j_v : vec_map_matrix[i])
					{
						if (abs(j_v.second) > ERROR_FOR_SM)
						{
							tri.push_back(Triplet<double>(i, j_v.first, j_v.second));
						}
					}
				}
				sparse_matrix.setFromTriplets(tri.begin(), tri.end());
			}
			/************************************************************************/
			/* 从 vector< map<int, double> > 转换到 Matrix矩阵                      */
			/************************************************************************/
			void TransVecMap2Matrix(const vector< map<int, double> >& vec_map_matrix, MatrixXd& matrix)
			{
				int nCount = (int)(vec_map_matrix.size());
				matrix.resize(nCount, nCount);
				matrix.setZero();
				for (int i = 0; i < (int)(vec_map_matrix.size()); ++i)
				{
					for (const pair<int, double>& j_v : vec_map_matrix[i])
					{
						matrix(i, j_v.first) = j_v.second;
					}
				}
			}
		} //	end of namespace PERIDYNAMIC
	} // end of namespace SAE
} // end of namespace DLUT

#endif