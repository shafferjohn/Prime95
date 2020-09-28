// WorkerDlg.h : header file
//
// Copyright 1995-2020 Mersenne Research, Inc.  All rights reserved
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
	CButton	c_workgroup;
	CStatic	c_threadnum_text;
	CComboBox c_threadnum;
	CStatic	c_work_pref_text;
	CComboBox c_work_pref;
	int	m_work_pref[MAX_NUM_WORKER_THREADS];
	CStatic	c_numcpus_text;
	CEdit	c_numcpus;
	int	m_numcpus[MAX_NUM_WORKER_THREADS];
	CButton	c_cert_work;
	int	m_cert_work;
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
	afx_msg void OnCbnKillfocusThreadNum();
	afx_msg void OnCbnKillfocusWorkType();
	afx_msg void OnEnKillfocusNumCpus();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

