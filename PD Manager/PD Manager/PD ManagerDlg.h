
// PD ManagerDlg.h : 头文件
//

#pragma once


// CPDManagerDlg 对话框
class CPDManagerDlg : public CDialogEx
{
// 构造
public:
	CPDManagerDlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_PDMANAGER_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonPdInformation();
	afx_msg void OnBnClickedRadioImplicit();
	afx_msg void OnBnClickedRadioExplicit();	
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedRadioConstantHorizon();
	afx_msg void OnBnClickedRadioConstantM();
	afx_msg void OnBnClickedRadioModal();
	CString m_fp_information;
	double 	m_ratio_of_delta_meshsize;
	double	m_time_step;
	int		m_plot_frames;
	int		m_radio_analysis_type;
	int		m_iterator_nums;
	int 	m_load_step;
	int		m_radio_horizon_type;
	double	m_d_const_horizon;
};
