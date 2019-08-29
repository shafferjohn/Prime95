// WorkerDlg.h : header file
//
// Copyright 1995-2017 Mersenne Research, Inc.  All rights reserved
//

/////////////////////////////////////////////////////////////////////////////
// CWorkerDlg dialog

class CWorkerDlg : public CDialog
{
// Construction
public:
	CWorkerDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CWorkerDlg)
	enum { IDD = IDD_WORKER_THREADS };
	CStatic	c_num_thread_text;
	CEdit	c_num_thread;
	UINT	m_num_thread;
	UINT	m_priority;
	CButton	c_workgroup;
	CStatic	c_threadnum_text;
	CComboBox c_threadnum;
	CStatic	c_work_pref_text;
	CComboBox c_work_pref;
	int	m_work_pref[MAX_NUM_WORKER_THREADS];
	CStatic	c_numcpus_text;
	CEdit	c_numcpus;
	int	m_numcpus[MAX_NUM_WORKER_THREADS];
	CButton c_hyper_tf;
	BOOL	m_hyper_tf;
	CButton c_hyper_ll;
	BOOL	m_hyper_ll;
	CStatic c_warn1_text;
	CStatic c_warn2_text;
	CStatic c_warn3_text;
	//}}AFX_DATA

private:
	int	thread_num;
public:
	int AreAllTheSame (int *);
	void InitComboBoxText (void);

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CWorkerDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CWorkerDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnEnKillfocusNumThread();
	afx_msg void OnEnSetfocusPriority();
	afx_msg void OnEnKillfocusPriority();
	afx_msg void OnCbnKillfocusThreadNum();
	afx_msg void OnCbnKillfocusWorkType();
	afx_msg void OnEnKillfocusNumCpus();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

