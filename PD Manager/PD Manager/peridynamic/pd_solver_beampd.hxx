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
	
					initializeMatParasPD();
				}
				/************************************************************************/
				/* ��ʽ��� F=Kd ƽ�ⷽ��                                               */
				/************************************************************************/
				bool				ImplicitSolve(int current_step, int total_step, int& ITERATOR_NUMS, int delta_step = 1)
				{
					TPdModel& pdModel = *m_pPdModel;

					ITERATOR_NUMS = 0;
					bool CONVERGENCED = FALSE;

					const set<int>& eids = pdModel.PdMeshCore().GetElementIdsByAll();
					const set<int>& nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					const int nCount = (int)(nids.size());

					/************************************************************************/
					/* �ڵ���Ϣ��ʼ��												         */
					/************************************************************************/
					parallel_for_each(nids.begin(), nids.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							node.IteratorDisplacement().setZero();
							node.IncrementalDisplacement().setZero();
						});
					/************************************************************************/
					/* ������٤�ɽ�Ԫ�ڲ����ֵ㴦����Ϣ��ʼ��                             */
					/************************************************************************/
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							for (int is = 0; is < IP_COUNT_2D; ++is)
							{
								element.IP(is).IteratorDisplacement().setZero();
								element.IP(is).IncrementalDisplacement().setZero();
							}
						});
					/************************************************************************/
					/* P                                                           */
					/************************************************************************/
					vector<double> P;
					P.clear();
					P.resize(DOF * nCount, 0);

					for (int nid : nids)
					{
						const TPdNode& node = pdModel.PdMeshCore().Node(nid);
						const vector<TLoadNodePoint>& LPs = node.LoadNodePoint();
						for (const TLoadNodePoint& lp : LPs)
						{
							double value = 0;
							int curid = lp.Lcid();

							if (pdModel.CurveExist(curid))
							{
								TCurve& curve = pdModel.Curve(curid);
								value = (curve.GetValueByX(current_step) - curve.GetValueByX(0)) * lp.Sf();
							}
							else
							{
								double cur_index = (double)(current_step) / (double)(total_step);
								value = lp.Sf() * cur_index;
							}

							int row = DOF * nid + (lp.Dof() - 1);
							P[row] = value;
						}
					}
					/************************************************************************/
					/* �������                                          */
					/************************************************************************/
					while (ITERATOR_NUMS < MAX_INTERATOR_NUMS && !CONVERGENCED)
					{
						/************************************************************************/
						/* R                                          */
						/************************************************************************/
						vector<double> R;
						R.clear();
						R.resize(DOF * nCount, 0);
						for (int ni : nids)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(ni);
							for (int loop_dim = 0; loop_dim < DOF; ++loop_dim)
							{
								R[ni * DOF + loop_dim] = node.InnerForce()[loop_dim];
							}
						}
						/************************************************************************/
						/* F=P-R                                          */
						/************************************************************************/
						vector<double> FORCE;
						FORCE.clear();
						FORCE.resize(DOF * nCount, 0);
						for (int loop = 0; loop < DOF * nCount; ++loop)
						{
							FORCE[loop] = P[loop] - R[loop];
						}
						/************************************************************************/
						/* �����ܸվ���                                                         */
						/************************************************************************/
						m_GK.clear();
						m_GK.resize(DOF * nCount);
						genGlobalStiffnessFEM();
						genGlobalStiffnessPD();

						/************************************************************************/
						/* ʩ��Լ������                                                         */
						/************************************************************************/
						for (int nid : nids)
						{
							const TPdNode& node = pdModel.PdMeshCore().Node(nid);
							const vector<TBoundarySpcNode>& BSNs = node.BoundarySpcNode();
							for (const TBoundarySpcNode& tbsn : BSNs)
							{
								int curSeri = DOF * nid + (tbsn.Dof() - 1);

								for (const pair<int, double>& j_v : m_GK[curSeri])
								{
									if (j_v.first != curSeri)
									{
										m_GK[j_v.first].erase(curSeri);
									}
								}

								m_GK[curSeri].clear();
								m_GK[curSeri][curSeri] = 1;

								FORCE[curSeri] = 0;
							}
						}
						/************************************************************************/
						/* ʩ��ǿ��λ��                                                         */
						/************************************************************************/
						for (int nid : nids)
						{
							const TPdNode& node = pdModel.PdMeshCore().Node(nid);
							const vector<TBoundaryPrescribedMotion>& BPMs = node.BoundaryPreMotion();
							for (const TBoundaryPrescribedMotion& bpm : BPMs)
							{
								int curid = bpm.Lcid();
								TCurve& curve = pdModel.Curve(curid);
								// ֻ��λ�Ʊ߽��������д���
								if (bpm.Vda() == 2)
								{
									if (ITERATOR_NUMS == 0)
									{
										//	�������߽�������λ�Ƶļ���
										double b_current = curve.GetValueByX((double)current_step / (double)total_step) * bpm.Sf();
										double b_last = curve.GetValueByX((double)(current_step - delta_step) / (double)(total_step)) * bpm.Sf();
										double b = b_current - b_last;

										int k = DOF * nid + (bpm.Dof() - 1);
										double Krr = m_GK[k][k];
										m_GK[k][k] = Krr * 10E10;

										FORCE[k] = Krr * b * 10E10;
									}
									else
									{
										int curSeri = DOF * nid + (bpm.Dof() - 1);

										for (const pair<int, double>& j_v : m_GK[curSeri])
										{
											if (j_v.first != curSeri)
											{
												m_GK[j_v.first].erase(curSeri);
											}
										}

										m_GK[curSeri].clear();
										m_GK[curSeri][curSeri] = 1;

										FORCE[curSeri] = 0;
									}
								}
							}
						}

						/************************************************************************/
						/* ������Է�����
						*/
						/************************************************************************/
						SparseMatrix<double> sparse_matrix_GK;
						TransVecMap2SparseMatrix(m_GK, sparse_matrix_GK);

						vector<double> DISPLACEMENT;
						DISPLACEMENT.clear();
						DISPLACEMENT.resize(DOF * nCount, 0);

						umf_solver(sparse_matrix_GK, DISPLACEMENT, FORCE);
						/************************************************************************/
						/* ����λ���������							                            */
						/************************************************************************/
						for (int nid : nids)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							for (int loop_dim = 0; loop_dim < DOF; ++loop_dim)
							{
								node.IteratorDisplacement()[loop_dim] = DISPLACEMENT[DOF * nid + loop_dim];
								node.IncrementalDisplacement()[loop_dim] += DISPLACEMENT[DOF * nid + loop_dim];
							}
						}

						/************************************************************************/
						/* ����λ�ƽ�������ж��Ƿ��Ѿ�����                                     */
						/************************************************************************/
					//	CONVERGENCED = isConvergenced();
						CONVERGENCED = true;
						if (CONVERGENCED)
						{
							/************************************************************************/
							/* ��������֮�󣬸��½ڵ������λ����Ϣ������������ѭ��                 */
							/************************************************************************/
							updateInfoAfterConvergenceFEM();
							updateInfoAfterConvergencePD();

							//		RefreshFracture(); 

							break;
						}
						else
						{
							/************************************************************************/
							/* ������������ֻ����Ӧ��/Ӧ��/��Ԫ&�ڵ���������������һ�ε���          */
							/************************************************************************/
							updateStrainStressFEM();
							updateStrainStressPD();

							parallel_for_each(nids.begin(), nids.end(), [&](int nid)
								{
									TPdNode& node = pdModel.PdMeshCore().Node(nid);
									node.InnerForce().setZero();
								});
							updateInnerForceFEM();
							updateInnerForcePD();
						}

						ITERATOR_NUMS++;
					}

					return CONVERGENCED;
				}
				/************************************************************************/
				/* ��ʽ��� ���Ĳ�ַ�                                                  */
				/************************************************************************/
				void				ExplicitSolve(double time_interval)
				{
					TPdModel& pdModel = *m_pPdModel;
					calPdForces();

					const set<int>& nodeIds = pdModel.PdMeshCore().GetNodeIdsByAll();
					parallel_for_each(nodeIds.begin(), nodeIds.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);

							int nodePos = node.Id() * DOF;
							node.Acceleration().x() = (node.OuterForce().x() - node.InnerForce().x()) / m_GM[nodePos][nodePos];
							node.Acceleration().y() = (node.OuterForce().y() - node.InnerForce().y()) / m_GM[nodePos + 1][nodePos + 1];
							node.Acceleration().z() = (node.OuterForce().z() - node.InnerForce().z()) / m_GM[nodePos + 2][nodePos + 2];
							node.Acceleration().rx() = (node.OuterForce().rx() - node.InnerForce().rx()) / m_GM[nodePos + 3][nodePos + 3];
							node.Acceleration().ry() = (node.OuterForce().ry() - node.InnerForce().ry()) / m_GM[nodePos + 4][nodePos + 4];
							node.Acceleration().rz() = (node.OuterForce().rz() - node.InnerForce().rz()) / m_GM[nodePos + 5][nodePos + 5];

							node.Velocity() = node.Velocity() + node.Acceleration() * time_interval;
							node.IteratorDisplacement() = node.Velocity() * time_interval;
							node.Displacement() = node.Displacement() + node.IteratorDisplacement();
						});
				}
				/************************************************************************/
				/* ���¶���ʧЧ��Bond��Ϣ                                               */
				/************************************************************************/
				void				RefreshFracture()
				{
					TPdModel& pdModel = *m_pPdModel;
					
					updateStrainStressPD();
					const set<int>& eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							if (element_i.AnalysisElementType() == PD_ELEMENT)
							{
								const double& sed_criterion = element_i.CalParas().sed_criterion;
								const Eigen::Matrix4d& D_PD = element_i.CalParas().D_PD;
								MAP_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
								for (MAP_NJ_FAMILY_ELEMENT::iterator iter = familyElements.begin();
									iter != familyElements.end(); ++iter)
								{								
									TPdFamilyElement& family_elem = iter->second;
									TPdElement& element_j = pdModel.PdMeshCore().Element(family_elem.Id());
									//	��Ԫ�������ж����ж�
									if (element_i.Id() == element_j.Id())
									{
										continue;
									}
									bool needUpdateSK = false;

									vector<TPdBond>& bonds = family_elem.Bonds();
									int nCountBondFractured = 0;
									for (TPdBond& bond : bonds)
									{
										if (bond.IsValid() && (bond.MicroPotential() > sed_criterion))
										{
											bond.MakeFailure();
											needUpdateSK = true;
										}
									}
									//	�����familyElement��bond������ʧЧ������Ҫ���¶�Ӧ��SK
									if (needUpdateSK)
									{
										genSingleStiffnessPD_Element(ei, family_elem);
									}
								}
							}
						}); 
				}
			private:
				/************************************************************************/
				/* Begin of PD															*/
				/************************************************************************/
				void				initializeMatParasPD()
				{
					TPdModel& pdModel = *m_pPdModel;

					double start = clock();

					for (const TPart& part : pdModel.Parts())
					{
						const TMaterial& material = pdModel.Material(part.MaterialId());
						const TSection& section = pdModel.Section(part.SectionId());
						string stype = section.Type();

						/************************************************************************/
						/*  ��ȡ���ϲ���                                                        */
						/************************************************************************/					
						double E, PR, Rho, Gc;
						if (material.Name() == "MAT_ELASTIC")
						{
							E= material.GetMatValue("E");
							PR = material.GetMatValue("PR");
							Rho = material.GetMatValue("Rho");
							Gc = material.GetMatValue("STRESS_TENSILE");
						}
						double t = section.GetSectionValue("THICKNESS");

						/************************************************************************/
						/* ����΢ģ��                                                           */
						/************************************************************************/
						const set<int> eids = part.GetElementIds();
						parallel_for_each(eids.begin(), eids.end(), [&](int ei) {
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
											
							double horizon = 0.0;
							if (DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON)
							{
								horizon = DLUT::SAE::PERIDYNAMIC::HORIZON;
							}
							else
							{
								horizon = element_i.SideLength() * DLUT::SAE::PERIDYNAMIC::RATIO_OF_HORIZON_MESHSIZE;
							}
							Eigen::Matrix4d& D_PD = element_i.CalParas().D_PD;
							double& sed_criterion = element_i.CalParas().sed_criterion;

							double D0 = E * pow(t, 3) / (12.0 * (1 - pow(PR, 2)));
							D_PD(0, 0) = 6 * E / (PI * pow(horizon, 3) * (1 - PR) * t);
							D_PD(1, 1) = 6 * D0 * (1 + PR) / (PI * pow(horizon, 3) * pow(t, 2));
							D_PD(2, 2) = E * (1 - 3 * PR) / (6 * PI * horizon * (1 - pow(PR, 2)) * t);
							D_PD(3, 3) = 6 * D0 * (1 - 3 * PR) / (PI * pow(horizon, 3) * pow(t, 2));

							sed_criterion = (3 * Gc) / (2 * t * pow(horizon, 3));
						});	
					}

					double total_time = (clock() - start) / 1000;
					cout << "initializeMatParasPD():\t\t" << total_time << endl;

					genSingleStiffnessPD();
					genGlobalMassPD();					
				}	
				void				genSingleStiffnessPD()
				{
					double start = clock();
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					/*parallel_*/for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							MAP_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
							for (MAP_NJ_FAMILY_ELEMENT::iterator iter = familyElements.begin();
								iter != familyElements.end(); ++iter)
							{
								TPdFamilyElement& family_elem = iter->second;
								genSingleStiffnessPD_Element(ei, family_elem);
							}
						});

					double total_time = (clock() - start) / 1000;
					cout << "genSingleStiffnessPD(): \t" << total_time << endl;
				}
				void				genSingleStiffnessPD_Element(int ei, TPdFamilyElement& family_elem)
				{
					TPdModel& pdModel = *m_pPdModel;
					TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
					Eigen::Matrix4d& D_PD = element_i.CalParas().D_PD;

					double thickness_i = element_i.Thickness();
					{
						int ej = family_elem.Id();
						const TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
						double thickness_j = element_j.Thickness();

						double ia = 0.5 * Distance_2pt<Vector3d>(element_i.Node(0).Coordinate(), element_i.Node(1).Coordinate());
						double ib = 0.5 * Distance_2pt<Vector3d>(element_i.Node(1).Coordinate(), element_i.Node(2).Coordinate());;

						double ja = 0.5 * Distance_2pt<Vector3d>(element_j.Node(0).Coordinate(), element_j.Node(1).Coordinate());
						double jb = 0.5 * Distance_2pt<Vector3d>(element_j.Node(1).Coordinate(), element_j.Node(2).Coordinate());;

						SingleStiffness& SK = family_elem.SK();
						SK.setZero();

						double volume_scale = family_elem.VolumeIndex();
						MatrixXd Te = T_elem(element_i, element_j);
						MatrixXd TE = T_ELEM(element_i, element_j);

						//	����ÿһ��bond�ĵ���
						const vector<TPdBond>& bonds = family_elem.Bonds();
						for (const TPdBond& bond : bonds)
						{
							if (bond.IsValid())
							{
								const TBondPoint& Xi = bond.Xi();
								const TBondPoint& Xj = bond.Xj();
								double L = bond.BondLength();

								double Ji = element_i.Jacobi(S_IP_2D[Xi.Index()], T_IP_2D[Xi.Index()]);
								double Jj = element_j.Jacobi(S_IP_2D[Xj.Index()], T_IP_2D[Xj.Index()]);
								
								// Bond ����Ϊ12*12�ľ���
								MatrixXd k;
								k.resize(12, 12);
								k.setZero();
								for (int js = 0; js < IP_COUNT_1D; ++js)
								{
									Eigen::MatrixXd BL = BL_Matrix_1D(L, S_IP_1D[js]);
									k += BL.transpose() * D_PD * BL * H_IP_1D[js] * L;
								}
								// Bond ����

								Eigen::MatrixXd Nij = N_IJ(S_IP_2D[Xi.Index()], T_IP_2D[Xi.Index()], ia, ib, S_IP_2D[Xj.Index()], T_IP_2D[Xj.Index()], ja, jb);
								Eigen::MatrixXd Tb = bond.T_b();

								//	0.5��ʾbond�����������һ��
								SK += 0.5 * volume_scale * Ji * thickness_i * Jj * thickness_j * TE.transpose() * Nij.transpose() * Te * Tb.transpose() * k * Tb * Te.transpose() * Nij * TE;
							}
						}
					}
				}

				void				genGlobalStiffnessPD()
				{
					TPdModel& pdModel = *m_pPdModel;
					const set<int>& nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					const set<int>& eids = pdModel.PdMeshCore().GetElementIdsByAll();

					/************************************************************************/
					/* ��װ�ܸվ���                                                         */
					/************************************************************************/
					for (int ei : eids)
					{
						const TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
						//	����ԪI�͵�ԪJ��Ӧ�Ľڵ����һ��nids������Ϊ8
						const vector<int>& nids_i = element_i.NodeIds();

						const MAP_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
						for (const PAIR_NJ_FAMILY_ELEMENT& familyElemInfo : familyElements)
						{
							int ej = familyElemInfo.first;
							const TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
							const vector<int>& nids_j = element_j.NodeIds();

							vector<int> nids(nids_i);
							for (int nj : nids_j)
							{
								nids.push_back(nj);
							}

							const SingleStiffness& SK = familyElemInfo.second.SK();
							int LenOfSK = 48;

							for (int row = 0; row < LenOfSK; ++row)
							{
								for (int col = 0; col < LenOfSK; ++col)
								{
									int Row = nids[row / DOF] * DOF + row % DOF;
									int Col = nids[col / DOF] * DOF + col % DOF;
									m_GK[Row][Col] += SK(row, col);
								}
							}
						}
					}
				}
				void				updateStrainStressPD()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					//	���·�����٤�ɽ�Ԫ���ֵ㴦��λ��������Ϣ
					parallel_for_each(eids.begin(), eids.end(), [&](int eid)
						{
							TPdElement& element = pdModel.PdMeshCore().Element(eid);
							int nid1 = element.NodeId(0);
							int nid2 = element.NodeId(1);
							int nid3 = element.NodeId(2);
							int nid4 = element.NodeId(3);
							const double a = element.SideLength();
							const double b = element.SideLength();

							const TPdNode& node1 = pdModel.PdMeshCore().Node(nid1);
							const TPdNode& node2 = pdModel.PdMeshCore().Node(nid2);
							const TPdNode& node3 = pdModel.PdMeshCore().Node(nid3);
							const TPdNode& node4 = pdModel.PdMeshCore().Node(nid4);

							MatrixXd delta_u_global;
							delta_u_global.resize(24, 1);
							delta_u_global.block(0, 0, 6, 1) = node1.IteratorDisplacement();
							delta_u_global.block(6, 0, 6, 1) = node2.IteratorDisplacement();
							delta_u_global.block(12, 0, 6, 1) = node3.IteratorDisplacement();
							delta_u_global.block(18, 0, 6, 1) = node4.IteratorDisplacement();

							MatrixXd TE;
							TE.resize(24, 24);
							TE.setZero();
							TE.block(0, 0, 3, 3) = element.LocalCoorSystem();
							TE.block(3, 3, 3, 3) = element.LocalCoorSystem();
							TE.block(6, 6, 6, 6) = TE.block(0, 0, 6, 6);
							TE.block(12, 12, 6, 6) = TE.block(0, 0, 6, 6);
							TE.block(18, 18, 6, 6) = TE.block(0, 0, 6, 6);
						
							MatrixXd delta_u_local = TE * delta_u_global;

							MatrixXd TN;
							TN.resize(6, 6);
							TN.setZero();
							TN.block(0, 0, 3, 3) = element.LocalCoorSystem();
							TN.block(3, 3, 3, 3) = element.LocalCoorSystem();

							double L = element.SideLength();
							//	���µ�Ԫ���ֵ㴦��λ��������ȫ������ϵ��
							//	��Ҫע�⣺Ҫ�Ե�Ԫ���Ĵ��Ļ��ֵ������Ϣ���£���Ϊ��ԪI����������
							for (int is = 0; is <= IP_COUNT_2D; ++is)
							{
								Eigen::MatrixXd N = M_SF_RECTANGLE_SHELL(S_IP_2D[is], T_IP_2D[is], a, b);
								//	ͨ��TN���ֲ�����ϵ�µĽڵ����������������ȫ������ϵ��
								element.IP(is).IteratorDisplacement() = TN.transpose() * N * delta_u_local;
								element.IP(is).IncrementalDisplacement() += element.IP(is).IteratorDisplacement();
							}
						});

					//	����bond���ֵ㴦��Ӧ���Ӧ��
					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
							//	�Դ�FEM����������bond�������Ϣ
							if (element_i.AnalysisElementType() != FEM_ELEMENT)
							{
								double alpha_i = element_i.Alpha();
								const Eigen::Matrix4d& D_PD = element_i.CalParas().D_PD;
								MAP_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
								for (MAP_NJ_FAMILY_ELEMENT::iterator iter = familyElements.begin();
									iter != familyElements.end(); ++iter)
								{
									int nj = (*iter).first;
									TPdFamilyElement& family_elem = (*iter).second;
									const TPdElement& element_j = pdModel.PdMeshCore().Element(nj);
									double alpha_j = element_j.Alpha();
									double alpha = (alpha_i + alpha_j) / 2.0;
									double volume_scale = family_elem.VolumeIndex();
					
									vector<TPdBond>& bonds = family_elem.Bonds();
									for (TPdBond& bond : bonds)
									{
										const TBondPoint& Xi = bond.Xi();
										const TBondPoint& Xj = bond.Xj();

										double L = bond.BondLength();

										MatrixXd delta_u_global;
										delta_u_global.resize(12, 1);
										delta_u_global.block(0, 0, 6, 1) = Xi.IteratorDisplacement();
										delta_u_global.block(6, 0, 6, 1) = Xj.IteratorDisplacement();

										MatrixXd T;
										T.resize(12, 12);
										T.setZero();
										T.block(0, 0, 3, 3) = bond.LocalCoorSystem();
										T.block(3, 3, 3, 3) = bond.LocalCoorSystem();
										T.block(6, 6, 3, 3) = bond.LocalCoorSystem();
										T.block(9, 9, 3, 3) = bond.LocalCoorSystem();

										MatrixXd delta_u_local = T * delta_u_global;
 										for (int js = 0; js < IP_COUNT_1D; ++js)
										{
											Eigen::MatrixXd BL = BL_Matrix_1D(L, S_IP_1D[js]);
											Eigen::MatrixXd BN_star = BN_star_Matrix_1D(L, S_IP_1D[js], delta_u_local);

											TStrain_Bond delta_strain_local = (BL + BN_star) * delta_u_local;
											TStress_Bond delta_stress_local = D_PD * delta_strain_local;

											bond.MicroPotential() += ((bond.IP(js).StressCurrent() + delta_stress_local * 0.5).transpose() * delta_strain_local)(0, 0) * H_IP_1D[js] * L;
											
											bond.IP(js).StrainCurrent() += delta_strain_local;
											bond.IP(js).StressCurrent() += delta_stress_local;
										}
									}
								}
							}
						});
				}
				void				updateInnerForcePD()
				{

				}
				void				updateInfoAfterConvergencePD()
				{

				}
				void	 			calPdForces()
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

					/************************************************************************/
					/* ���㵥ԪI�뵥ԪJ֮���������  48*1 ����                              */
					/************************************************************************/
					const set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					parallel_for_each(eids.begin(), eids.end(), [&](int ei)
						{
							TPdElement& element_i = pdModel.PdMeshCore().Element(ei);

							MAP_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
							for (MAP_NJ_FAMILY_ELEMENT::iterator iter = familyElements.begin();
								iter != familyElements.end(); ++iter)
							{
								int ej = iter->first;
								TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
							
								//	����ԪI�͵�ԪJ��Ӧ�Ľڵ����һ��nids������Ϊ8							
								vector<int> nids = element_i.NodeIds();								
								for (int nj : element_j.NodeIds())
								{
									nids.push_back(nj);
								}
								const SingleStiffness& SK = iter->second.SK();

								VectorXd DIS;
								DIS.resize(48);
								DIS.setZero();
								int loop_Dis = 0;
								for (int nid : nids)
								{
									DIS.block(6*loop_Dis++, 0, 6, 1) = pdModel.PdMeshCore().Node(nid).Displacement();
								}
								iter->second.ForceOfBond() = SK * DIS;													
							}
						});

					/************************************************************************/
					/* �����е��������ۼӵ���Ԫ�ڵ���ȥ                                     */
					/************************************************************************/
					for (int ei : eids)
					{
						TPdElement& element_i = pdModel.PdMeshCore().Element(ei);
						MAP_NJ_FAMILY_ELEMENT& familyElements = element_i.FamilyElements();
						for (const PAIR_NJ_FAMILY_ELEMENT& familyElemInfo : familyElements)
						{
							int ej = familyElemInfo.first;
							const TPdElement& element_j = pdModel.PdMeshCore().Element(ej);
							vector<int> nids = element_i.NodeIds();
							const vector<int>& nids_j = element_j.NodeIds();							
							for (int nj : nids_j)
							{
								nids.push_back(nj);
							}
							const VectorXd& FORCE = familyElemInfo.second.ForceOfBond();
							
							int loop_force = 0;
							for (int ni : nids)
							{
								pdModel.PdMeshCore().Node(ni).InnerForce() += FORCE.block(6 * loop_force++, 0, 6, 1);
							}
						}
					}
				}
				void				genGlobalMassPD()
				{
					double start = clock();
					TPdModel& pdModel = *m_pPdModel;

					int nCount = pdModel.PdMeshCore().NodeCount();
					m_GM.clear();
					m_GM.resize(DOF * nCount);

					for (const TPart& part : pdModel.Parts())
					{
						const set<int> eids = part.GetElementIds();
						const TMaterial& material = pdModel.Material(part.MaterialId());
						double Rho = material.GetMatValue("Rho");

						/************************************************************************/
						/* ��װ��������                                                         */
						/************************************************************************/
						for (int ei : eids)
						{
							const TPdElement& element = pdModel.PdMeshCore().Element(ei);
							//	����ԪI�͵�ԪJ��Ӧ�Ľڵ����һ��nids������Ϊ8
							const vector<int>& nids = element.NodeIds();
							double h = element.Thickness();
							double A = element.Area();
							for (int ni : nids)
							{
								int CurPos = ni * DOF;
								m_GM[CurPos][CurPos] += Rho * A * h / 4.0;
								m_GM[CurPos + 1][CurPos + 1] += Rho * A * h / 4.0;
								m_GM[CurPos + 2][CurPos + 2] += Rho * A * h / 4.0;

								m_GM[CurPos + 3][CurPos + 3] += Rho * A * pow(h, 3) / 24.0;
								m_GM[CurPos + 4][CurPos + 4] += Rho * A * pow(h, 3) / 24.0;
								m_GM[CurPos + 5][CurPos + 5] += Rho * A * pow(h, 3) / 12.0;
							}
						}
					}
					double total_time = (clock() - start) / 1000;
					cout << "genGlobalMassPD():\t\t" << total_time << endl;
				}		
			private:
				Eigen::MatrixXd N_IJ(double is, double it, double ia, double ib, double js, double jt, double ja, double jb)
				{
					Eigen::MatrixXd Res;
					Res.resize(12, 48);
					Res.setZero();

					Res.block(0, 0, 6, 24) = M_SF_RECTANGLE_SHELL(is, it, ia, ib);
					Res.block(6, 24, 6, 24) = M_SF_RECTANGLE_SHELL(js, jt, ja, jb);

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
				Eigen::MatrixXd SK_LOCAL(double cax, double cby, double cbz, double ctor, double L)
				{
					MatrixXd k;
					k.resize(12, 12);
					k.setZero();

					k(0, 0) = cax / L;
					k(0, 6) = -k(0, 0);

					k(1, 1) = 12.0 * cbz / pow(L, 3);
					k(1, 5) = 6.0 * cbz / pow(L, 2);
					k(1, 7) = -k(1, 1);
					k(1, 11) = k(1, 5);

					k(2, 2) = 12.0 * cby / pow(L, 3);
					k(2, 4) = -6.0 * cby / pow(L, 2);
					k(2, 8) = -k(2, 2);
					k(2, 10) = k(2, 4);

					k(3, 3) = ctor / L;
					k(3, 9) = -k(3, 3);

					k(4, 2) = k(2, 4);
					k(4, 4) = 4.0 * cby / L;
					k(4, 8) = -k(4, 2);
					k(4, 10) = 2.0 * cby / L;

					k(5, 1) = k(1, 5);
					k(5, 5) = 4.0 * cbz / L;
					k(5, 7) = -k(5, 1);
					k(5, 11) = 2.0 * cbz / L;

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
					k(10, 10) = 4.0 * cby / L;

					k(11, 1) = k(1, 11);
					k(11, 5) = k(5, 11);
					k(11, 7) = k(7, 11);
					k(11, 11) = 4.0 * cbz / L;

					return k;
				}

			private:
				/************************************************************************/
				/* Begin of FEM															*/
				/************************************************************************/
				void				genSingleStiffnessFEM()
				{

				}
				void				genGlobalStiffnessFEM()
				{

				}
				void				updateStrainStressFEM()
				{

				}
				void				updateInnerForceFEM()
				{

				}
				void				updateInfoAfterConvergenceFEM()
				{
					TPdModel& pdModel = *m_pPdModel;
					set<int> eids = pdModel.PdMeshCore().GetElementIdsByAll();
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();

					//	����FEM�ڵ���Ϣ
					parallel_for_each(nids.begin(), nids.end(), [&](int nid)
						{
							TPdNode& node = pdModel.PdMeshCore().Node(nid);
							/************************************************************************/
							/* ��ǰ��������λ���ۼӵ�ȫ�ֽڵ���λ��                                 */
							/* �ڵ㵱ǰ����λ�Ƹ���Ϊ��һ��λ��+��ǰ��������λ��                    */
							/************************************************************************/
							node.Displacement() += node.IncrementalDisplacement();
							node.Coordinate() += node.IncrementalDisplacement().block(0, 0, 3, 1);
						});
				}
				/************************************************************************/
				/* End of FEM															*/
				/************************************************************************/
								
				/************************************************************************/
				/* Begin of PD&FEM														*/
				/************************************************************************/
				Eigen::MatrixXd		G_Matrix_1D(double L, double chi)
				{
					Eigen::MatrixXd Res;
					Res.resize(2, 12);
					Res.setZero();

					Res(0, 2) = 6 * (chi * chi - chi) / L;
					Res(0, 4) = -(3 * chi * chi - 4 * chi + 1);
					Res(0, 8) = -Res(0, 2);
					Res(0, 10) = -(3 * chi * chi - 2 * chi);

					Res(1, 1) = 6 * (chi * chi - chi) / L;
					Res(1, 5) = (3 * chi * chi - 4 * chi + 1);
					Res(1, 7) = -Res(0, 2);
					Res(1, 11) = (3 * chi * chi - 2 * chi);

					return Res;
				}
				Eigen::MatrixXd		BL_Matrix_1D(double L, double chi)
				{
					Eigen::MatrixXd Res;
					Res.resize(4, 12);
					Res.setZero();

					Res(0, 0) = -1.0 / L;
					Res(0, 6) = 1.0 / L;

					Res(1, 2) = (12 * chi - 6) / (L * L);
					Res(1, 4) = -(6 * chi - 4) / L;
					Res(1, 8) = -Res(1, 2);
					Res(1, 10) = -(6 * chi - 2) / L;

					Res(2, 1) = (12 * chi - 6) / (L * L);
					Res(2, 5) = (6 * chi - 4) / L;
					Res(2, 7) = -Res(2, 1);
					Res(2, 11) = (6 * chi - 2) / L;

					Res(3, 3) = -1.0 / L;
					Res(3, 9) = 1.0 / L;

					return Res;
				}
				Eigen::MatrixXd		BN_star_Matrix_1D(double L, double chi, const MatrixXd& delta_u)
				{
					Eigen::MatrixXd Res;
					Res.resize(4, 12);
					Res.setZero();

					//Eigen::MatrixXd G = G_Matrix_1D(L, chi);

					//Eigen::MatrixXd delta_A;
					//delta_A.resize(4, 2);
					//delta_A.setZero();
					//delta_A.block(0, 0, 1, 2) = (G * delta_u).transpose();

					//Res = 0.5 * delta_A * G;

					return Res;
				}

				bool				isConvergenced()
				{
					bool CONVERGENCED = true;

					TPdModel& pdModel = *m_pPdModel;
					set<int> nids = pdModel.PdMeshCore().GetNodeIdsByAll();
					double MaxOfIterU = -1E20;
					double MaxOfIncrU = -1E20;

					double MaxOfIterDisp = -1E20;

					for (int nid : nids)
					{
						TPdNode& node = pdModel.PdMeshCore().Node(nid);
						const vector<TBoundaryPrescribedMotion>& BPMs = node.BoundaryPreMotion();
						//	ʩ��ǿ��λ�ƵĽڵ㲻��Ϊƽ���жϵĵ�
						if (BPMs.size() == 0)
						{
							double IterDisp = Module(node.IteratorDisplacement());
							if (IterDisp > MaxOfIterDisp)
							{
								MaxOfIterDisp = IterDisp;
							}

							//for (int loop_dim = 0; loop_dim < DOF; ++loop_dim)
							//{
							//	double IterU = abs(node.IteratorDisplacement()[loop_dim]);
							//	double IncreU = abs(node.IncrementalDisplacement()[loop_dim]);
							//	if (IncreU > MaxOfIncrU)
							//	{
							//		MaxOfIncrU = IncreU;
							//		MaxOfIterU = IterU;
							//	}
							//}
						}
					}


					/*double ratio = MaxOfIterU / MaxOfIncrU;
					if (ratio > CONVERGENCE_FACTOR)*/
					if (MaxOfIterDisp > CONVERGENCE_FACTOR)
					{
						CONVERGENCED = false;
					}
					else
					{
						CONVERGENCED = true;
					}

					return CONVERGENCED;
				}
			private:			
				vector< map<int, double> >	m_GK;	//	Stiffness matrix
				vector< map<int, double> >  m_GM;	//	Mass matrix
			private:
				TPdModel*					m_pPdModel;
			};
		} // END OF PERIDYNAMIC
	} // END OF DLUT
} // END OF SAE

#endif