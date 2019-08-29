// StartDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CStartDlg dialog

class CStartDlg : public CDialog
{
// Construction
public:
	CStartDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CStartDlg)
	enum { IDD = IDD_WORKER_START };
	CStatic	c_thread_text;
	CEdit	c_thread;
	UINT	m_thread;
	BOOL	m_all_threads;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CStartDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CStartDlg)
		// NOTE: the ClassWizard will add member functions here
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

	afx_msg void OnBnClickedRadio1();
};
