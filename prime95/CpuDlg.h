// CpuDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CCpuDlg dialog

class CCpuDlg : public CDialog
{
// Construction
public:
	CCpuDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
	//{{AFX_DATA(CCpuDlg)
	enum { IDD = IDD_CPU };
	UINT	m_hours;
	CString	m_start_time;
	CString	m_end_time;
	UINT	m_day_memory;
	UINT	m_night_memory;
	CString	m_cpu_info;
	BOOL	m_memory_editable;
	CStatic	c_day_memory_text;
	CEdit	c_day_memory;
	CStatic	c_night_memory_text;
	CEdit	c_night_memory;
	CStatic	c_start_time_text;
	CEdit	c_start_time;
	CStatic	c_end_time_text;
	CEdit	c_end_time;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CCpuDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CCpuDlg)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};
