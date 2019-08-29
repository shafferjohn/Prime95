// StopDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CStopDlg dialog

class CStopDlg : public CDialog
{
// Construction
public:
	CStopDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CStopDlg)
	enum { IDD = IDD_WORKER_STOP };
	CStatic	c_thread_text;
	CEdit	c_thread;
	UINT	m_thread;
	BOOL	m_all_threads;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CStopDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CStopDlg)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

	afx_msg void OnBnClickedRadio1();
};
