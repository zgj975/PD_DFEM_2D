#ifndef DLUT_SAE_PERIDYNAMIC_PD_SOLVERS_BEAM_PD_HXX_20201225
#define DLUT_SAE_PERIDYNAMIC_PD_SOLVERS_BEAM_PD_HXX_20201225

#include "pd_database.hxx"

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			//	PD Implicit analysis 
			class TSolverBeamPD
			{
			public:
				TSolverBeamPD() {}
				~TSolverBeamPD() {}
			public:
				void				Attach(TPdModel& model)
				{
					m_pPdModel = &model;

					initializeMatParas();
				}
				/************************************************************************/
				/* ��ʽ��� F=Kd ƽ�ⷽ��                                               */
				/************************************************************************/
				void				ImplicitSolve(int current_step, int total_step, int delta_step)
				{
					TPdModel& pdModel = *m_pPdModel;

					/************************************************************************/
					/* �����ܸվ���                                                         */
					/************************************************************************/
					double start = clock();					
					genGlobalStiffnessPD();
					cout << "genGlobalStiffnessPD():\t\t" << (clock() - start) / 1000. << endl;

					/************************************************************************/
					/* ����������                                                           */
					/************************************************************************/
					start = clock();
					int nCount = pdModel.PdMeshCore().NodeCount();
			
					int dim = DOF;
					m_GF.clear();
					m_GF.resize(dim * nCount, 0);

					for (const TLoadNodePoint& lp : pdModel.LoadNodePoints())
					{
						/************************************************************************/
						/* �������ص�ֻ����������˽���������ƽ���ֲ���ÿ����������             */
						/************************************************************************/
						double incr_index = (double)(delta_step) / (double)(total_step);

						int nid = lp.Id();

						double value = 0;
						int curid = lp.Lcid();

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
					/* ʩ��Լ������                                                         */
					/************************************************************************/
					for (const TBoundarySpcNode& tbsn : pdModel.BoundarySpcNodes())
					{
						int nid = tbsn.Id();
						
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
					/* ʩ��ǿ��λ��                                                         */
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
							// ֻ��λ�Ʊ߽��������д���
							if (bpm.Vda() == 2)
							{
								//	�������߽�������λ�Ƶļ���
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

					cout << "Boundary Conditions:\t\t" << (clock() - start) / 1000. << endl;
					start = clock();
					
					/************************************************************************/
					/* ������Է�����						                                */
					/************************************************************************/
					SparseMatrix<double> sparse_matrix_GK;
					TransVecMap2SparseMatrix(m_GK, sparse_matrix_GK);
					
					cout << "TransVecMap2SparseMatrix():\t" << (clock() - start) / 1000. << endl;
					start = clock();

					m_GD.clear();
					m_GD.resize(dim * nCount, 0);
		//			cout << sparse_matrix_GK << endl;
					umf_solver(sparse_matrix_GK, m_GD, m_GF);
					/************************************************************************/
					/* ����λ�ƽ��							                                */
					/************************************************************************/
					start = clock();
					for (int nid = 0; nid < nCount; ++nid)
					{
						TPdNode& node = pdModel.PdMeshCore().Node(nid);

						node.Displacement().x() += m_GD[dim * nid];
						node.Displacement().y() += m_GD[dim * nid + 1];
						node.Displacement().z() += m_GD[dim * nid + 2];
						node.Displacement().rx() += m_GD[dim * nid + 3];
						node.Displacement().ry() += m_GD[dim * nid + 4];
						node.Displacement().rz() += m_GD[dim * nid + 5];
					}
	
					cout << "umf_solver:\t\t\t" << (clock() - start) / 1000. << endl;
				}
				/************************************************************************/
				/* ��ʽ��� ���Ĳ�ַ�                                                  */
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
					}
				}
				/************************************************************************/
				/* ���¶���ʧЧ��Bond��Ϣ                                               */
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
								/* ���¶Ե�Ԫi���е��յ�����                                           */
								/************************************************************************/
								if (updateSK)
								{
									genSingleStiffnessPD(element_i);
								}
							}
						});
					}
				}
			private:
				/************************************************************************/
				/* ��ʼ��PD������ļ������                                             */
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
						/*  ��ȡ���ϲ���                                                        */
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
						/* ����΢ģ��                                                           */
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

							double& cax = element_i.CalParas().cax;
							double& cby = element_i.CalParas().cby;
							double& cbz = element_i.CalParas().cbz;
							double& ctor = element_i.CalParas().ctor;
							double& csy = element_i.CalParas().csy;

							double ri = element_i.Radius();
							double D0 = E * pow(t, 3) / (12.0 * (1 - pow(PR, 2)));
							cax = 6 * E / (PI * (pow(horizon, 3) - pow(ri, 3)) * (1 - PR) * t);
							//	���Ǽ���Ӱ��
							double phi = 5 * (horizon - ri) / (2 * dx * sqrt(15 * (1 + PR)));
							double lamda = phi / (phi - atan(phi));
							cbz = lamda * E * (1 - 3 * PR) / (6 * PI * (horizon - ri) * (1 - pow(PR, 2)) * t);
												
							cby = 6 * D0 * (1 + PR) / (PI * (pow(horizon, 3) - pow(ri, 3)) * pow(t, 2));
							ctor = 6 * D0 * (1 - 3 * PR) / (PI * (pow(horizon, 3) - pow(ri, 3)) * pow(t, 2));

							double G = E / (2 * (1 + PR));
							double k = 6.0 / 5.0;
							csy = 6 * G / (k * PI * t * (pow(horizon, 3) - pow(ri, 3)));
						});
						/************************************************************************/
						/* ��������Bond�ĵ��վ���                                               */
						/************************************************************************/
						double start = clock();
						parallel_for_each(eids.begin(), eids.end(), [&](int ei)
							{
								TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
								genSingleStiffnessPD(element_i);
							});
						double total_time = (clock() - start) / 1000;
						cout << "genSingleStiffnessPD():\t\t" << total_time << endl;
					}
				}	
				/************************************************************************/
				/* �����ܸվ���		   		                                            */
				/************************************************************************/
				void				genGlobalStiffnessPD()
				{
					TPdModel& pdModel = *m_pPdModel;

					int nCount = pdModel.PdMeshCore().NodeCount();
					int dim = DOF;
					m_GK.clear();
					m_GK.resize(dim * nCount);

					for (const TPart& part : pdModel.Parts())
					{
						const set<int> eids = part.GetElementIds();
						const TSection& section = pdModel.Section(part.SectionId());
						string stype = section.Type();
							
						/************************************************************************/
						/* ��װ�ܸվ���                                                         */
						/************************************************************************/
						for (int ei : eids)
						{
							const TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							//	����ԪI�͵�ԪJ��Ӧ�Ľڵ����һ��nids������Ϊ8
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
								int LenOfSK = 48;
							
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
				/* ���ɵ�ԪI�뵥ԪJ֮��ĵ��վ���                                       */
				/************************************************************************/
				void				genSingleStiffnessPD(TPdElement& element_i)
				{
					TPdModel& pdModel = *m_pPdModel;
					TPart& part = pdModel.Part(element_i.PartId());
					TSection& section = pdModel.Section(part.SectionId());
					double PR = pdModel.Material(part.MaterialId()).GetMatValue("PR");
				
					const double cax = element_i.CalParas().cax;
					const double cby = element_i.CalParas().cby;
					const double cbz = element_i.CalParas().cbz;
					const double ctor = element_i.CalParas().ctor;
					const double csy = element_i.CalParas().csy;

					map<int, TPdBond>& familyBonds = element_i.FamilyElementBonds();
					for (map<int, TPdBond>::iterator iter = familyBonds.begin();
						iter != familyBonds.end(); ++iter)
					{
						int nj = (*iter).first;
						TPdBond& bond_ij = (*iter).second;
						const TPdElement& element_j = pdModel.PdMeshCore().Element(nj);

						SingleStiffness& SK = bond_ij.SK();
						SK.resize(48, 48);
						SK.setZero();

						double volume_scale = bond_ij.Volume() / (element_j.Area() * element_j.Thickness());

						double s[2], t[2];
						s[0] = t[0] = 1.0 / sqrt(3.0);
						s[1] = t[1] = -s[0];												

						MatrixXd Te = T_elem(element_i, element_j);
						MatrixXd TE = T_ELEM(element_i, element_j);

						for (int is = 0; is < 2; ++is)
						{
							for (int it = 0; it < 2; ++it)
							{
								TCoordinate Xi = element_i.CoordinateInElement(s[is], t[it]);
								double thickness_i = element_i.Thickness();
								double Ji = element_i.Jacobi(s[is], t[it]);
								double ia = 0.5 * Distance_2pt<Vector3d>(element_i.Node(0).CoordinateCurrent(), element_i.Node(1).CoordinateCurrent());
								double ib = 0.5 * Distance_2pt<Vector3d>(element_i.Node(1).CoordinateCurrent(), element_i.Node(2).CoordinateCurrent());;

								for (int js = 0; js < 2; ++js)
								{
									for (int jt = 0; jt < 2; ++jt)
									{
										TCoordinate Xj = element_j.CoordinateInElement(s[js], t[jt]);
										double thickness_j = element_j.Thickness();
										double Jj = element_j.Jacobi(s[js], t[jt]);
										
										double ja = 0.5 * Distance_2pt<Vector3d>(element_j.Node(0).CoordinateCurrent(), element_j.Node(1).CoordinateCurrent());
										double jb = 0.5 * Distance_2pt<Vector3d>(element_j.Node(1).CoordinateCurrent(), element_j.Node(2).CoordinateCurrent());;

										double L = Module<TCoordinate>(Xj - Xi);
										MatrixXd k = SK_LOCAL(cax, cby, cbz, ctor, csy, L, element_i.SideLength(), PR);
										Eigen::MatrixXd Nij = N_IJ(s[is], t[it],ia, ib, s[js], t[jt], ja, jb);

										//	����Ҫ���й�һ������
										Eigen::Vector3d vec_ij = Xj - Xi;
										Eigen::MatrixXd Tb = T_b(vec_ij);

										// 0.5��ʾbond�����������һ��
										SK += 0.5 * volume_scale * Ji * Jj * TE.transpose() * Nij.transpose() * Te * Tb.transpose() * k * Tb * Te.transpose() * Nij * TE * thickness_i * thickness_j;
									}
								}
							}
						}
					}
				}
				/************************************************************************/
				/* ����PD��֮������ܶ�	                                                */
				/************************************************************************/
				void				calPdForces()
				{
					TPdModel& pdModel = *m_pPdModel;

					/************************************************************************/
					/* ÿһ������Ҫ�����������ʼ��                                         */
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
							//	����ԪI�͵�ԪJ��Ӧ�Ľڵ����һ��nids������Ϊ8
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
								DIS.resize(48);
								DIS.setZero();
								int loop_Dis = 0;
								for (int ni : nids_i)
								{
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().x();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().y();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().z();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().rx();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().ry();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(ni).Displacement().rz();
								}
								for (int nj : nids_j)
								{
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().x();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().y();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().z();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().rx();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().ry();
									DIS[loop_Dis++] = pdModel.PdMeshCore().Node(nj).Displacement().rz();
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
								pdModel.PdMeshCore().Node(ni).InnerForce().z() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(ni).InnerForce().rx() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(ni).InnerForce().ry() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(ni).InnerForce().rz() += FORCE[loop_force++];
							}
							for (int nj : nids_j)
							{
								pdModel.PdMeshCore().Node(nj).InnerForce().x() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(nj).InnerForce().y() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(nj).InnerForce().z() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(nj).InnerForce().rx() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(nj).InnerForce().ry() += FORCE[loop_force++];
								pdModel.PdMeshCore().Node(nj).InnerForce().rz() += FORCE[loop_force++];
							}
						}
					}
				}
			private:
				Eigen::MatrixXd N_IJ(double is, double it, double ia, double ib, double js, double jt, double ja, double jb)
				{
					Eigen::MatrixXd Res;
					Res.resize(12, 48);
					Res.setZero();

					double pI[4] = { -1, 1, 1, -1 };
					double qI[4] = { -1, -1, 1, 1 };

					Res.block(0, 0, 6, 6) = N_Ipq(is, it, pI[0], qI[0], ia, ib);
					Res.block(0, 6, 6, 6) = N_Ipq(is, it, pI[1], qI[1], ia, ib);
					Res.block(0, 12, 6, 6) = N_Ipq(is, it, pI[2], qI[2], ia, ib);
					Res.block(0, 18, 6, 6) = N_Ipq(is, it, pI[3], qI[3], ia, ib);

					Res.block(6, 24, 6, 6) = N_Ipq(js, jt, pI[0], qI[0], ja, jb);
					Res.block(6, 30, 6, 6) = N_Ipq(js, jt, pI[1], qI[1], ja, jb);
					Res.block(6, 36, 6, 6) = N_Ipq(js, jt, pI[2], qI[2], ja, jb);
					Res.block(6, 42, 6, 6) = N_Ipq(js, jt, pI[3], qI[3], ja, jb);

					return Res;
				}
				Eigen::MatrixXd T_elem(const TPdElement& element_i, const TPdElement& element_j)
				{
					Eigen::MatrixXd Res;
					Res.resize(12, 12);
					Res.setZero();

					Matrix3d Local_Ei = element_i.LocalCoorSystem();
					Matrix3d Local_Ej = element_j.LocalCoorSystem();

					Res.block(0, 0, 3, 3) = Local_Ei;
					Res.block(3, 3, 3, 3) = Local_Ei;
					Res.block(6, 6, 3, 3) = Local_Ej;
					Res.block(9, 9, 3, 3) = Local_Ej;

					return Res;
				}
				Eigen::MatrixXd T_ELEM(const TPdElement& element_i, const TPdElement& element_j)
				{
					Eigen::MatrixXd Res;
					Res.resize(48, 48);
					Res.setZero();

					Matrix3d Local_Ei = element_i.LocalCoorSystem();
					Matrix3d Local_Ej = element_j.LocalCoorSystem();

					for (int loop = 0; loop < 8; ++loop)
					{
						Res.block(loop * 3, loop * 3, 3, 3) = Local_Ei;
						Res.block((loop + 8) * 3, (loop + 8) * 3, 3, 3) = Local_Ej;
					}

					return Res;
				}
				Eigen::Matrix3d LAMDA(const Vector3d& localVector)
				{
					Eigen::Matrix3d T;
					T.setZero();

					double idist = Module<Vector3d>(localVector);
					double Cx = localVector.x() / idist;
					double Cy = localVector.y() / idist;
					double Cz = localVector.z() / idist;

					if (Cz + ERR_VALUE > 1.0)
					{
						//	�ֲ�x����ȫ��Z��һ��
						T(0, 2) = 1;
						T(1, 1) = 1;
						T(2, 0) = -1;
					}
					else if (Cz - ERR_VALUE < -1.0)
					{
						//	�ֲ�x����ȫ��Z�᷽���෴
						T(0, 2) = -1;
						T(1, 1) = 1;
						T(2, 0) = 1;
					}
					else
					{
						double D = sqrt(Cx * Cx + Cy * Cy);
						T(0, 0) = Cx;
						T(0, 1) = Cy;
						T(0, 2) = Cz;

						T(1, 0) = -Cy / D;
						T(1, 1) = Cx / D;
						T(1, 2) = 0;

						T(2, 0) = -Cx * Cz / D;
						T(2, 1) = -Cy * Cz / D;
						T(2, 2) = D;
					}
					return T;
				}
				Eigen::MatrixXd T_b(const Vector3d& localVector)
				{
					Eigen::MatrixXd Res;
					Res.resize(12, 12);
					Res.setZero();
					Eigen::Matrix3d lamda = LAMDA(localVector);
					Res.block(0, 0, 3, 3) = lamda;
					Res.block(3, 3, 3, 3) = lamda;
					Res.block(6, 6, 3, 3) = lamda;
					Res.block(9, 9, 3, 3) = lamda;

					return Res;
				}
				Eigen::MatrixXd N_Ipq(double p, double q, double pI, double qI, double a, double b)
				{
					Eigen::MatrixXd Res;
					Res.resize(6, 6);
					Res.setZero();

					Res(0, 0) = (1 + p * pI) * (1 + q * qI) / 4.0;
					Res(1, 1) = Res(0, 0);

					Res(2, 2) = -((1 + pI * p) * (1 + qI * q) * (p * p + q * q - pI * p - qI * q - 2)) / 8.0;
					Res(2, 3) = b * qI * ((1 + pI * p) * (1 + qI * q) * (1 + qI * q) * (qI * q - 1))/ 8.0;
					Res(2, 4) = -a * pI * ((1 + pI * p) * (1 + pI * p) * (1 + qI * q) * (pI * p - 1)) / 8.0;

					Res(3, 2) = (-qI * (pI * p + 1) * (3 * q * q + p * p - pI * p - 3)) / (b * 8.0);
					Res(3, 3) = (b * (3 * qI * q - 1) * (1 + pI * p) * (1 + qI * q)) / (b * 8.0);
					Res(3, 4) = (-a * pI * qI * (1 + pI * p) * (1 + pI * p) * (pI * p - 1)) / (b * 8.0);
		
					Res(4, 2) = (-pI * (qI * q + 1) * (3 * p * p + q * q - qI * q -3))/ (-a * 8.0);
					Res(4, 3) = (b * pI * qI * (qI * q + 1) * (qI * q + 1) * (qI * q - 1)) / (-a * 8.0);
					Res(4, 4) = (-a * (3 * pI * p - 1) * (pI * p + 1) * (qI * q + 1)) / (-a * 8.0);

					Res(5, 5) = Res(0, 0);

					return Res;
				}
				Eigen::MatrixXd SK_LOCAL(double cax, double cby, double cbz, double ctor, double csy, double L, double dx, double PR)
				{
					//double By = 12 * cby / (pow(L, 2) * csy);
					//double Bz = 12 * (1 + PR) * pow(dx, 2) / (5 * pow(L, 2));
					double By = 0;
					double Bz = 0;

					MatrixXd k;
					k.resize(12, 12);
					k.setZero();

					k(0, 0) = cax / L;
					k(0, 6) = -k(0, 0);

					k(1, 1) = 12.0 * cbz / pow(L, 3) / (1 + Bz);
					k(1, 5) = 6.0 * cbz / pow(L, 2) / (1 + Bz);
					k(1, 7) = -k(1, 1);
					k(1, 11) = k(1, 5);

					k(2, 2) = 12.0 * cby / pow(L, 3) / (1 + By);
					k(2, 4) = -6.0 * cby / pow(L, 2) / (1 + By);
					k(2, 8) = -k(2, 2);
					k(2, 10) = k(2, 4);

					k(3, 3) = ctor / L;
					k(3, 9) = -k(3, 3);

					k(4, 2) = k(2, 4);
					k(4, 4) = (4.0 + By) * cby / L / (1 + By);
					k(4, 8) = -k(4, 2);
					k(4, 10) = (2.0 - By) * cby / L / (1 + By);

					k(5, 1) = k(1, 5);
					k(5, 5) = (4.0 + Bz) * cbz / L / (1 + Bz);
					k(5, 7) = -k(5, 1);
					k(5, 11) = (2.0 - Bz) * cbz / L / (1 + Bz);

					k(6, 0) = k(0, 6);
					k(6, 6) = -k(6, 0);

					k(7, 1) = k(1, 7);
					k(7, 5) = k(5, 7);
					k(7, 7) = -k(7, 1);
					k(7, 11) = k(7, 5);

					k(8, 2) = k(2, 8);
					k(8, 4) = k(4, 8);
					k(8, 8) = -k(8, 2);
					k(8, 10) = k(8, 4);

					k(9, 3) = k(3, 9);
					k(9, 9) = -k(9, 3);

					k(10, 2) = k(2, 10);
					k(10, 4) = k(4, 10);
					k(10, 8) = k(8, 10);
					k(10, 10) = (4.0 + By) * cby / L / (1 + By);

					k(11, 1) = k(1, 11);
					k(11, 5) = k(5, 11);
					k(11, 7) = k(7, 11);
					k(11, 11) = (4.0 + Bz) * cbz / L / (1 + Bz);

					return k;
				}
			private:
				// F=K*D
				vector<double>				m_GF;	//	Force vector
				vector<double>				m_GD;	//	Displacement vector
				vector< map<int, double> >	m_GK;	//	Stiffness matrix
			private:
				TPdModel*					m_pPdModel;
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif