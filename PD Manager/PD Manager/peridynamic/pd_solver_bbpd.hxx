/*==============================================================================

Copyright 2017 Dalian University of Technology .
All rights reserved

================================================================================
-- Please append file description informations here --
The Solver for 2D blank / shell
================================================================================
Date            Name                    Description of Change
2020 / 04 / 21		Zheng Guojun			Create
2020 / 04 / 21		Zheng Guojun			Update RefreshFracture
$HISTORY$
================================================================================
*/

#ifndef DLUT_SAE_PERIDYNAMIC_PD_SOLVERS_BBPD_HXX_20200421
#define DLUT_SAE_PERIDYNAMIC_PD_SOLVERS_BBPD_HXX_20200421

#include "pd_database.hxx"

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			//	PD Implicit analysis 
			class TSolverBBPD
			{
			public:
				TSolverBBPD() {}
				~TSolverBBPD() {}
			public:
				void				Attach(TPdModel& model)
				{
					m_pPdModel = &model;

					initializeMatParas();
				}
				/************************************************************************/
				/* 隐式求解 F=Kd 平衡方程                                               */
				/************************************************************************/
				void				ImplicitSolve(int current_step, int total_step, int delta_step)
				{
					TPdModel& pdModel = *m_pPdModel;

					/************************************************************************/
					/* 生成总刚矩阵                                                         */
					/************************************************************************/
					double start = clock();
					cout << "Begin to Initialize Material Parameters of Peridynamic model..." << endl;
					
					genGlobalStiffness();
					cout << "Stiffness Matrix:\t\t\t\t" << (clock() - start) / 1000. << endl;

					/************************************************************************/
					/* 生成力向量                                                           */
					/************************************************************************/
					start = clock();
					int nCount = pdModel.PdMeshCore().NodeCount();
			
					int dim = 2;
					m_GF.clear();
					m_GF.resize(dim * nCount, 0);

					for (const TLoadNodePoint& lp : pdModel.LoadNodePoints())
					{
						/************************************************************************/
						/* 静力加载的只有总力，因此将力的增量平均分布到每个增量步中             */
						/************************************************************************/
						double incr_index = (double)(delta_step) / (double)(total_step);

						int nid = lp.Id();

						double value = 0;
						int curid = lp.Lcid();
						if (lp.Dof() == 3 || lp.Dof() == 4 || lp.Dof() == 5 || lp.Dof() == 6)
						{
							continue;
						}

						if (pdModel.CurveExist(curid))
						{
							TCurve& curve = pdModel.Curve(curid);
							value = (curve.GetValueByX(current_step) - curve.GetValueByX(current_step - delta_step)) * lp.Sf();
						}
						else
						{
							value = lp.Sf() * incr_index;
						}

						int row = dim * nid + (lp.Dof() - 1);
						m_GF[row] = value;
					}

					/************************************************************************/
					/* 施加约束条件                                                         */
					/************************************************************************/
					for (const TBoundarySpcNode& tbsn : pdModel.BoundarySpcNodes())
					{
						int nid = tbsn.Id();
						if (tbsn.Dof() == 3 || tbsn.Dof() == 4 || tbsn.Dof() ==5 || tbsn.Dof() == 6)
						{
							continue;
						}
						int curSeri = dim * nid + (tbsn.Dof() - 1);

						for (const pair<int, double>& j_v : m_GK[curSeri])
						{
							if (j_v.first != curSeri)
							{
								m_GK[j_v.first].erase(curSeri);
							}
						}

						m_GK[curSeri].clear();
						m_GK[curSeri][curSeri] = 1;

						m_GF[curSeri] = 0;
					}

					/************************************************************************/
					/* 施加强制位移                                                         */
					/************************************************************************/
					for (const TBoundaryPrescribedMotion& bpm : pdModel.BoundaryPrescribedMotions())
					{
						string type = bpm.MotionType();
						set<int>	nids;
						nids.clear();
						if (type == "RIGID")
						{
							int partId = bpm.Id();
							TPart& part = pdModel.Part(partId);
							nids = part.GetElementIds();
						}
						else if (type == "NODE")
						{
							nids.insert(bpm.Id());
						}
						for (int nid : nids)
						{
							int curid = bpm.Lcid();
							TCurve& curve = pdModel.Curve(curid);
							// 只对位移边界条件进行处理
							if (bpm.Vda() == 2)
							{
								//	按照曲线进行增量位移的计算
								double b_current = curve.GetValueByX((double)current_step / (double)total_step) * bpm.Sf();
								double b_last = curve.GetValueByX((double)(current_step - delta_step) / (double)(total_step)) * bpm.Sf();
								double b = b_current - b_last;

								if (bpm.Dof() == 3 || bpm.Dof() == 4 || bpm.Dof() == 5 || bpm.Dof() == 6)
								{
									continue;
								}

								int k = dim * nid + (bpm.Dof() - 1);
								double Krr = m_GK[k][k];
								m_GK[k][k] = Krr * 10E10;

								m_GF[k] = Krr * b * 10E10;
							}
						}
					}

					cout << "set boundary:\t\t\t\t\t" << (clock() - start) / 1000. << endl;
					
					/************************************************************************/
					/* 求解线性方程组						                                */
					/************************************************************************/
					SparseMatrix<double> sparse_matrix_GK;
					TransVecMap2SparseMatrix(m_GK, sparse_matrix_GK);
					
					cout << "transform to Sparse Matrix:\t\t\t" << (clock() - start) / 1000. << endl;
					start = clock();

					m_GD.clear();
					m_GD.resize(dim * nCount, 0);

					umf_solver(sparse_matrix_GK, m_GD, m_GF);
					/************************************************************************/
					/* 更新位移结果							                                */
					/************************************************************************/
					start = clock();
					for (int nid = 0; nid < nCount; ++nid)
					{
						TPdNode& node = pdModel.PdMeshCore().Node(nid);

						node.Displacement().x() += m_GD[dim * nid];
						node.Displacement().y() += m_GD[dim * nid + 1];
					}
	
					cout << "Update strain and strss:\t\t\t" << (clock() - start) / 1000. << endl;
				}
				/************************************************************************/
				/* 显式求解 中心差分法                                                  */
				/************************************************************************/
				void				ExplicitSolve(double time_interval)
				{
					TPdModel& pdModel = *m_pPdModel;
					calPdForces();

					for (const TPart& part : pdModel.Parts())
					{
						const TMaterial& material = pdModel.Material(part.MaterialId());
						if (material.Name() == "MAT_RIGID")
							continue;

						const TSection& section = pdModel.Section(part.SectionId());
						double Rho = material.GetMatValue("Rho");
						double thickness = section.GetSectionValue("THICKNESS");

						const set<int>& eleIds = part.GetElementIds();
						set<int> nodeIds;
						nodeIds.clear();
						for (int eid : eleIds)
						{
							const TPdElement& element = pdModel.PdMeshCore().Element(eid);
							for (int nid : element.NodeIds())
							{
								nodeIds.insert(nid);
							}
						}

						parallel_for_each(nodeIds.begin(), nodeIds.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							double nVolume = 0;
							set<int> adjEleIds = node.AdjElementIds();
							for (int adjEle : adjEleIds)
							{
								const TPdElement& element = pdModel.PdMeshCore().Element(adjEle);
								nVolume += 0.25 * element.Area() * element.Thickness();
							}
							node.Acceleration() = (node.OuterForce() - node.InnerForce()) / (Rho * nVolume);
							node.Velocity() = node.Velocity() + node.Acceleration() * time_interval;
							node.Displacement() = node.Displacement() + node.Velocity() * time_interval;
						});

						/*parallel_for_each(nodeIds.begin(), nodeIds.end(), [&](int nid)
							{
								TPdNode& node = pdModel.PdMeshCore().Node(nid);
								double nVolume = 0;
								set<int> adjEleIds = node.AdjElementIds();
								for (int adjEle : adjEleIds)
								{
									const TPdElement& element = pdModel.PdMeshCore().Element(adjEle);
									nVolume += 0.25 * element.Area() * element.Thickness();
								}

								node.Acceleration() = (node.OuterForce() - node.InnerForce()) / (Rho * nVolume);


							});*/
					}
				}
				/************************************************************************/
				/* 更新断裂失效的Bond信息                                               */
				/************************************************************************/
				void				RefreshFracture()
				{
					TPdModel& pdModel = *m_pPdModel;
					for (const TPart& part : pdModel.Parts())
					{
						const TMaterial& material = pdModel.Material(part.MaterialId());
						double E = material.GetMatValue("E");
						const set<int> eids = part.GetElementIds();
	
						double s[2], t[2];
						s[0] = t[0] = 1.0 / sqrt(3.0);
						s[1] = t[1] = -s[0];

						parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							double s0 = element_i.CalParas().s0;
							bool updateSK = false;

							if (element_i.CalParas().b_facture)
							{
								map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
								for (map<int, TPdBond>::iterator iter = familyBonds.begin();
									iter != familyBonds.end();)
								{
									int ej = (*iter).first;
									TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
									if (element_j.CalParas().b_facture)
									{
										//const TCoordinate& Xi = element_i.CoordinateInElement(0, 0);
										//const TDisplacement& Di = element_i.DisplaceInElement(0, 0);

										//const TCoordinate& Xj = element_j.CoordinateInElement(0, 0);
										//const TDisplacement& Dj = element_j.DisplaceInElement(0, 0);

										//double idist = Distance_2pt<Vector3d>(Xi, Xj);
										//double nlength = Distance_2pt<Vector3d>(Xi + Di.block(0, 0, 3, 1), Xj + Dj.block(0, 0, 3, 1));
										//double s = abs((nlength - idist) / idist);
										//if (s > s0)
										//{
										//	++iter;
										//	element_i.DeleteFamilyElement(element_j.Id());
										////	element_j.DeleteFamilyElement(element_i.Id());
										//	updateSK = true;
										//}
										//else
										//{
										//	++iter;
										//}
										vector<double> countForS;
										countForS.clear();
										for (int is = 0; is < 2; ++is)
										{
											for (int it = 0; it < 2; ++it)
											{
												const TCoordinate& Xi = element_i.CoordinateInElement(s[is], t[it]);
												const TDisplacement& Di = element_i.DisplaceInElement(s[is], t[it]);

												for (int js = 0; js < 2; ++js)
												{
													for (int jt = 0; jt < 2; ++jt)
													{
														const TCoordinate& Xj = element_j.CoordinateInElement(s[js], t[jt]);
														const TDisplacement& Dj = element_j.DisplaceInElement(s[js], t[jt]);

														double idist = Distance_2pt<Vector3d>(Xi, Xj);
														double nlength = Distance_2pt<Vector3d>(Xi + Di.block(0, 0, 3, 1), Xj + Dj.block(0, 0, 3, 1));
														double s = abs((nlength - idist) / idist);
														if (s > s0)
														{
															countForS.push_back(s);
														}
													}
												}
											}
										}

										if (countForS.size() >= 16)
										{
											++iter;
											element_i.DeleteFamilyElement(element_j.Id());
											updateSK = true;
										}
										else
										{
											++iter;
										}
									}
									else
									{
										++iter;
									}
								}
								/************************************************************************/
								/* 重新对单元i进行单刚的生成                                           */
								/************************************************************************/
								if (updateSK)
								{
									genSingleStiffness(element_i);
								}
							}
						});
					}
				}
			private:
				/************************************************************************/
				/* 初始化PD点所需的计算参数                                             */
				/************************************************************************/
				void				initializeMatParas()
				{
					TPdModel& pdModel = *m_pPdModel;
						
					for (const TPart& part : pdModel.Parts())
					{
						const TMaterial& material = pdModel.Material(part.MaterialId());
						const TSection& section = pdModel.Section(part.SectionId());
						string stype = section.Type();

						/************************************************************************/
						/*  获取材料参数                                                        */
						/************************************************************************/					
						double E, PR, Rho, G0;
						if (material.Name() == "MAT_ELASTIC")
						{
							E= material.GetMatValue("E");
							PR = material.GetMatValue("PR");
							Rho = material.GetMatValue("Rho");
							G0 = material.GetMatValue("STRESS_TENSILE");
						}
						double t = section.GetSectionValue("THICKNESS");

						/************************************************************************/
						/* 计算微模量                                                           */
						/************************************************************************/
						const set<int> eids = part.GetElementIds();
						parallel_for_each(eids.begin(), eids.end(), [&](int ei) {
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);

							const double dx = element_i.SideLength();
							double horizon = 0.0;
							if (DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON)
							{
								horizon = DLUT::SAE::PERIDYNAMIC::HORIZON;
							}
							else
							{
								horizon = element_i.SideLength() * DLUT::SAE::PERIDYNAMIC::RATIO_OF_HORIZON_MESHSIZE;
							}

							double& density = element_i.CalParas().density;
							double& c = element_i.CalParas().c;
							double& s0 = element_i.CalParas().s0;

							double ri = element_i.Radius();
							c = 9 * E / (PI * (pow(horizon, 3) - pow(ri, 3)));
						//	c = 9 * E / (PI * (pow(horizon, 3)));

							s0 = sqrt((4 * PI * G0) / (9 * E * horizon));
						});

						/************************************************************************/
						/* 计算所有Bond的单刚矩阵                                               */
						/************************************************************************/
						parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							genSingleStiffness(element_i);
						});
					}
				}	
				/************************************************************************/
				/* 生成总刚矩阵		   		                                            */
				/************************************************************************/
				void				genGlobalStiffness()
				{
					TPdModel& pdModel = *m_pPdModel;

					int nCount = pdModel.PdMeshCore().NodeCount();
					int dim = 2;
					m_GK.clear();
					m_GK.resize(dim * nCount);

					for (const TPart& part : pdModel.Parts())
					{
						const set<int> eids = part.GetElementIds();
						const TSection& section = pdModel.Section(part.SectionId());
						string stype = section.Type();
							
						/************************************************************************/
						/* 组装总刚矩阵                                                         */
						/************************************************************************/
						for (int ei : eids)
						{
							const TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							//	将单元I和单元J对应的节点放入一个nids，长度为8
							const vector<int>& nids_i = element_i.NodeIds();

							const map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
							for (const pair<int, TPdBond>& bondInfo : familyBonds)
							{
								int ej = bondInfo.first;
								const TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
								const vector<int>& nids_j = element_j.NodeIds();

								vector<int> nids(nids_i);								
								for (int nj : nids_j)
								{
									nids.push_back(nj);
								}

								const SingleStiffness& SK = bondInfo.second.SK();
								int LenOfSK = 16;
							
								for (int row = 0; row < LenOfSK; ++row)
								{
									for (int col = 0; col < LenOfSK; ++col)
									{
										int Row = nids[row / dim] * dim + row % dim;
										int Col = nids[col / dim] * dim + col % dim;
										m_GK[Row][Col] += SK(row, col);
									}
								}
							}
						}
					}
				}
				/************************************************************************/
				/* 生成单元I与单元J之间的单刚矩阵                                       */
				/************************************************************************/
				void				genSingleStiffness(TPdElement& element_i)
				{
					TPdModel& pdModel = *m_pPdModel;
					TPart& part = pdModel.Part(element_i.PartId());
					TSection& section = pdModel.Section(part.SectionId());
					double thickness = section.GetSectionValue("THICKNESS");

					map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
					for (map<int, TPdBond>::iterator iter = familyBonds.begin();
						iter != familyBonds.end(); ++iter)
					{
						int nj = (*iter).first;
						TPdBond& bond_ij = (*iter).second;
						const TPdElement& element_j = pdModel.PdMeshCore().Element(nj);

						SingleStiffness& SK = bond_ij.SK();
						SK.resize(16, 16);
						SK.setZero();

						const double c = element_i.CalParas().c;
						double volume_scale = bond_ij.Volume() / (element_j.Area() * element_j.Thickness());

						double s[2], t[2];
						s[0] = t[0] = 1.0 / sqrt(3.0);
						s[1] = t[1] = -s[0];

						MatrixXd k;
						k.resize(2, 2);
						k.setZero();

						k(0, 0) = 1;
						k(0, 1) = -1;
						k(1, 0) = -1;
						k(1, 1) = 1;

						MatrixXd Te = Telem(element_i, element_j);
						MatrixXd TE = TELEM(element_i, element_j);

						for (int is = 0; is < 2; ++is)
						{
							for (int it = 0; it < 2; ++it)
							{
								TCoordinate Xi = element_i.CoordinateInElement(s[is], t[it]);
								double Ji = element_i.Jacobi(s[is], t[it]);

								for (int js = 0; js < 2; ++js)
								{
									for (int jt = 0; jt < 2; ++jt)
									{
										Eigen::MatrixXd Nij = N(s[is], t[it], s[js], t[jt]);

										TCoordinate Xj = element_j.CoordinateInElement(s[js], t[jt]);
										double Jj = element_j.Jacobi(s[js], t[jt]);

										double idist = Module<TCoordinate>(Xj - Xi);

										//	方向要进行归一化处理
										Eigen::Vector2d vec_ij = (Xj - Xi).block(0, 0, 2, 1) / idist;

										Eigen::MatrixXd Tb;
										Tb.resize(2, 4);
										Tb.setZero();

										Tb(0, 0) = Tb(1, 2) = vec_ij.x();
										Tb(0, 1) = Tb(1, 3) = vec_ij.y();
										// 0.5表示bond正反都会计算一次
										SK += 0.5 * volume_scale * c / idist * Ji * Jj * TE * Nij.transpose() * Te.transpose() * Tb.transpose() * k * Tb * Te * Nij * TE.transpose();
									}
								}
							}
						}
					}
				}
				/************************************************************************/
				/* 计算PD点之间的力密度	                                                */
				/************************************************************************/
				void				calPdForces()
				{
					TPdModel& pdModel = *m_pPdModel;

					/************************************************************************/
					/* 每一步都需要将内力置零初始化                                         */
					/************************************************************************/
					const set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					parallel_for_each(nids.begin(), nids.end(), [&](int ni)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(ni);
							node.InnerForce().setZero();
						});

					const set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							//	将单元I和单元J对应的节点放入一个nids，长度为8
							vector<int> nids_i = element_i.NodeIds();

							map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
							for (map<int, TPdBond>::iterator iter = familyBonds.begin();
								iter != familyBonds.end(); ++iter)
							{
								int ej = iter->first;
								TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
								vector<int> nids_j = element_j.NodeIds();
								const SingleStiffness& SK = iter->second.SK();

								VectorXd DIS;
								DIS.resize(16);
								DIS.setZero();
								int loop_Dis = 0;
								for (int ni : nids_i)
								{
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().x();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().y();
								}
								for (int nj : nids_j)
								{
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().x();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().y();
								}

								iter->second.ForceOfBond() = SK * DIS;
							}
						});

					for (int ei : eids)
					{
						TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
						vector<int> nids_i = element_i.NodeIds();
						map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
						for (const pair<int, TPdBond>& bondInfo : familyBonds)
						{
							int ej = bondInfo.first;
							TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
							const VectorXd& FORCE = bondInfo.second.ForceOfBond();
							vector<int> nids_j = element_j.NodeIds();

							int loop_force = 0;
							for (int ni : nids_i)
							{
								pdModel.PdMeshCore().Node(ni).InnerForce().x() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(ni).InnerForce().y() += FORCE[loop_force++];
							}
							for (int nj : nids_j)
							{
								pdModel.PdMeshCore().Node(nj).InnerForce().x() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(nj).InnerForce().y() += FORCE[loop_force++];
							}
						}
					}
				}
			private:
				Eigen::MatrixXd N(double is, double it, double js, double jt)
				{
					Eigen::MatrixXd Res;
					Res.resize(4, 16);
					Res.setZero();
					double Ni1 = (1 - is) * (1 - it) / 4;
					double Ni2 = (1 + is) * (1 - it) / 4;
					double Ni3 = (1 + is) * (1 + it) / 4;
					double Ni4 = (1 - is) * (1 + it) / 4;
					Res(0, 0) = Res(1, 1) = Ni1;
					Res(0, 2) = Res(1, 3) = Ni2;
					Res(0, 4) = Res(1, 5) = Ni3;
					Res(0, 6) = Res(1, 7) = Ni4;

					double Nj1 = (1 - js) * (1 - jt) / 4;
					double Nj2 = (1 + js) * (1 - jt) / 4;
					double Nj3 = (1 + js) * (1 + jt) / 4;
					double Nj4 = (1 - js) * (1 + jt) / 4;
					Res(2, 8) = Res(3, 9) = Nj1;
					Res(2, 10) = Res(3, 11) = Nj2;
					Res(2, 12) = Res(3, 13) = Nj3;
					Res(2, 14) = Res(3, 15) = Nj4;

					return Res;
				}
				Eigen::MatrixXd Telem(const TPdElement& element_i, const TPdElement& element_j)
				{
					TPdModel& pdModel = *m_pPdModel;

					Eigen::MatrixXd Res;
					Res.resize(4, 4);
					Res.setZero();

					Matrix3d Local_Ei = element_i.LocalCoorSystem();
					Matrix3d Local_Ej = element_j.LocalCoorSystem();

					Res.block(0, 0, 2, 2) = Local_Ei.block(0, 0, 2, 2);
					Res.block(2, 2, 2, 2) = Local_Ej.block(0, 0, 2, 2);

					return Res;
				}
				Eigen::MatrixXd TELEM(const TPdElement& element_i, const TPdElement& element_j)
				{
					TPdModel& pdModel = *m_pPdModel;
					Eigen::MatrixXd Res;
					Res.resize(16, 16);
					Res.setZero();

					Matrix3d Local_Ei = element_i.LocalCoorSystem();
					Matrix3d Local_Ej = element_j.LocalCoorSystem();

					Res.block(0, 0, 2, 2) = Local_Ei.block(0, 0, 2, 2);
					Res.block(2, 2, 2, 2) = Local_Ei.block(0, 0, 2, 2);
					Res.block(4, 4, 2, 2) = Local_Ei.block(0, 0, 2, 2);
					Res.block(6, 6, 2, 2) = Local_Ei.block(0, 0, 2, 2);

					Res.block(8, 8, 2, 2) = Local_Ej.block(0, 0, 2, 2);
					Res.block(10, 10, 2, 2) = Local_Ej.block(0, 0, 2, 2);
					Res.block(12, 12, 2, 2) = Local_Ej.block(0, 0, 2, 2);
					Res.block(14, 14, 2, 2) = Local_Ej.block(0, 0, 2, 2);

					return Res;
				}
	
			private:
				// F=K*D
				vector<double>				m_GF;	//	Force vector
				vector<double>				m_GD;	//	Displacement vector
				vector< map<int, double> >	m_GK;	//	Stiffness matrix
				vector< map<int, double> >	m_GM;	//	Mass matrix
			private:
				TPdModel*					m_pPdModel;
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif