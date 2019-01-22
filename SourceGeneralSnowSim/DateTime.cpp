// $Id: date_time.C 442 2006-03-07 07:58:51Z bjg $
#include "DateTime.h"
#include <stdlib.h>

//#ifndef LINUX
//DateTime NoDateTime(1753, 1, 1, 0, 0, 0);     // Equal NULL value //erikt
//#endif
// 86400 = 24*60*60

extern void skipDelim(istream&);

const char* Month[12] =
{
    "Januar", "Februar", "Mars", "April", "Mai", "Juni", "Juli",
    "August", "September", "Oktober", "November", "Desember"
};

const char* Month_eng[12] =
{
    "January", "February", "March", "April", "May", "June", "July",
    "August", "September", "October", "November", "December"
};

const char* Month_short[12] =
{
    "Jan", "Feb", "Mar", "Apr", "Mai", "Jun", "Jul",
    "Aug", "Sep", "Okt", "Nov", "Des"
};

const char* Month_eng_short[12] =
{
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
    "Aug", "Sep", "Oct", "Nov", "Dec"
};

int legalDateTime(yearType y, monthType m, dayType d,
                  short h, short min, short s)
{
    return (legalDate(y, m, d) && legalTime(h, min, s)) ?
           1 : 0;
}

//DateTime::DateTime(char *datetime)

DateTime::DateTime() : Date(), CTime()
{
    // NOP
}

DateTime::DateTime(const char *date, const char* time) : Date(date),
    CTime(time)
{
    // NOP
}

DateTime::DateTime(yearType yy, monthType mm, dayType dd,
                   hourType hh, minuteType mn, secondType ss) :
    Date(yy, mm, dd), CTime(hh, mn, ss)
{
    // NOP
}

DateTime::DateTime(const Date &d) : Date(d), CTime((hourType)0, 0, 0)
{
    // NOP
}

DateTime::DateTime(const DateTime& d) : Date(d.julnum), CTime(d.tsec)
{
    // NOP
}

// 199401011200
// 19940101120000
DateTime::DateTime(const char *datetime)
{
    char temp[5];
    temp[4] = '\0';
    char temp1[3];
    temp1[2] = '\0';

    yearType  y = atoi(strncpy(temp, datetime, 4));
    monthType m = atoi(strncpy(temp1, datetime + 4, 2));
    dayType   d = atoi(strncpy(temp1, datetime + 6, 2));
    Date::julnum = date2jday(d, m, y);

    short h, mi, sec = 0;
    h = atoi(strncpy(temp1, datetime + 8, 2));
    mi = atoi(strncpy(temp1, datetime + 10, 2));
    if (strlen(datetime) > 12)
    {
        sec = atoi(strncpy(temp1, datetime + 12, 2));
    }
    CTime::tsec = time2totSec(h, mi, sec);

}

char *nextdigits(char *str, int num)
{
    int i = 0, onlydigits = 1;
    int strmax = strlen(str);

    while (((i < num && onlydigits) || !isdigit(str[i])) && i < strmax)
    {
        i++;
        if (!isdigit(str[i]))
        {
            onlydigits = 0;
        }
    }

    if (!isdigit(str[i]))
    {
        return NULL;
    }
    return str + i;
}

DateTime::DateTime(char *datetime, int form)
{
    char temp[5];
    temp[4] = '\0';
    char temp1[3];
    temp1[2] = '\0';
    char *ptr1, *ptr2, *ptr3, *ptr4, *ptr5;
    yearType  y;
    monthType m;
    dayType d;
    short h = 0, mi = 0, sec = 0;

    if (!datetime)
    {
        return;
    }

    switch (form)
    {
    case 0:
        ptr1 = strstr(datetime, "/") + 1;
        if (!ptr1)
        {
            y = 0;
            m = 2;
            d = 30;
            return;
        }
        ptr2 = strstr(ptr1, "/") + 1;
        if (!ptr2)
        {
            y = 0;
            m = 2;
            d = 30;
            return;
        }
        sscanf(datetime, "%hu", &m);
        sscanf(ptr1, "%hu", &d);
        sscanf(ptr2, "%hu", &y);
        ptr3 = strstr(ptr2, " ") + 1;
        if (ptr3)
        {
            sscanf(ptr3, "%hd", &h);
            ptr4 = strstr(ptr3, ":") + 1;
            if (ptr4)
            {
                sscanf(ptr4, "%hd", &mi);
                ptr5 = strstr(ptr4, ":") + 1;
                if (ptr5)
                {
                    sscanf(ptr5, "%hd", &sec);
                }
            }
        }
        break;
    case 1:
        d = atoi(strncpy(temp1, datetime, 2));
        ptr1 = nextdigits(datetime, 2);
        if (!ptr1)
        {
            y = 0;
            m = 2;
            d = 30;
            return;
        }
        m = atoi(strncpy(temp1, ptr1, 2));
        ptr2 = nextdigits(ptr1, 2);
        if (!ptr2)
        {
            y = 0;
            m = 2;
            d = 30;
            return;
        }
        y = atoi(strncpy(temp, ptr2, 4));
        ptr3 = nextdigits(ptr2, 4);
        if (!ptr3)
        {
            break;
        }
        h = atoi(strncpy(temp1, ptr3, 2));
        ptr4 = nextdigits(ptr3, 2);
        if (!ptr4)
        {
            break;
        }
        mi = atoi(strncpy(temp1, ptr4, 2));
        ptr5 = nextdigits(ptr4, 2);
        if (!ptr5)
        {
            break;
        }
        sec = atoi(strncpy(temp1, ptr5, 2));
        /*
        y = atoi(strncpy(temp, datetime+6, 4));
        m = atoi(strncpy(temp1, datetime+3, 2));
        d = atoi(strncpy(temp1, datetime, 2));
        h = atoi(strncpy(temp1, datetime+11, 2));
        mi = atoi(strncpy(temp1, datetime+14, 2));
        sec = 0;
        */
        break;
    case 2:
        y = atoi(strncpy(temp, datetime, 4));
        m = atoi(strncpy(temp1, datetime + 4, 2));
        d = atoi(strncpy(temp1, datetime + 6, 2));
        h = atoi(strncpy(temp1, datetime + 8, 2));
        mi = atoi(strncpy(temp1, datetime + 10, 2));
        sec = 0;
        break;
    case 3:
        y = atoi(strncpy(temp, datetime, 4));
        m = atoi(strncpy(temp1, datetime + 5, 2));
        d = atoi(strncpy(temp1, datetime + 8, 2));
        h = atoi(strncpy(temp1, datetime + 11, 2));
        mi = atoi(strncpy(temp1, datetime + 14, 2));
        sec = 0;
        break;
    case 4:
        strncpy(temp, datetime + 6, 4);
        temp[4] = '\0';
        y = atoi(temp);
        strncpy(temp1, datetime + 3, 2);
        temp1[2] = '\0';
        m = atoi(temp1);
        strncpy(temp1, datetime, 2);
        temp1[2] = '\0';
        d = atoi(temp1);
        h = 0;
        mi = 0;
        sec = 0;
        break;
    }

    Date::julnum = date2jday(d, m, y);
    CTime::tsec = time2totSec(h, mi, sec);
}




void DateTime::now()
{
    Date::now();
    CTime::now();
}

int DateTime::legal() const
{
    return (Date::legal() && CTime::legal()) ? 1 : 0;
}

istream &operator>>(istream& is, DateTime& dt)
{
    char date[9], time[7], temp[5], temp1[3];
    date[8] = time[6] = temp[4] = temp1[2] = '\0';
    skipDelim(is);
    is.getline(date, 9, '/');
    skipDelim(is);
    is.getline(time, 7);

    dt.julnum = dt.date2jday((dayType)atoi(strncpy(temp1, date + 6, 2)),
                             (monthType)atoi(strncpy(temp1, date + 4, 2)),
                             (yearType)atoi(strncpy(temp, date, 4)));

    dt.tsec = dt.time2totSec(atoi(strncpy(temp1, time, 2)),
                             atoi(strncpy(temp1, time + 2, 2)),
                             atoi(strncpy(temp1, time + 4, 2)));
    return is;
}

DateTime& DateTime::operator=(const DateTime &dt2)
{
    julnum = dt2.julnum;
    tsec = dt2.tsec;
    return *this;
}

DateTime& DateTime::operator=(const Date &dt2)
{
    julnum = dt2.julnum;
    return *this;
}

DateTime operator+ (int sec, const DateTime &dt)
{
    if (sec < 0)
    {
        return operator-(dt, abs(sec));
    }
    else
        return DateTime(dt.julnum + (dt.tsec + sec) / 86400,
                        (dt.tsec + sec) % 86400);
}

DateTime operator+ (const DateTime &dt, int sec)
{
    return ::operator+(sec, dt);
}

void operator+= (DateTime &dt, int sec)
{
    if (sec < 0)
    {
        dt -= abs(sec);
    }
    else
    {
        dt.julnum += ((dt.tsec + sec) / (86400));
        dt.tsec = ((dt.tsec + sec) % (86400));
    }
}

DateTime DateTime::Floor(int precision_min)
{
    DateTime dt2(*this);
    int min = dt2.tsec / 60;

    min -= min % precision_min;

    dt2.tsec = min * 60;
    return dt2;
}

DateTime DateTime::Ceiling(int precision_min)
{
    DateTime dt2(*this);
    int min = dt2.tsec / 60;

    if (min % precision_min)
    {
        min -= min % precision_min;
        min += precision_min;

        if (min == 24 * 60)
        {
            min = 0;
            dt2 += 24 * 60 * 60;
        }

        dt2.tsec = min * 60;
    }

    return dt2;
}

unsigned long operator- (const DateTime &dt1, const DateTime &dt2)
{
    // this func. returns the absolute difference between the two times
    // in seconds

    unsigned long result;
    if (dt1 >= dt2)
    {
        if (dt1.tsec < dt2.tsec)
        {
            result = dt1.tsec + 86400 - dt2.tsec;
            result += (dt1.julnum - 1 - dt2.julnum) * 86400;
        }
        else
        {
            result = dt1.tsec - dt2.tsec;
            result += (dt1.julnum - dt2.julnum) * 86400;
        }
    }
    else
    {
        if (dt2.tsec < dt1.tsec)
        {
            result = dt2.tsec + 86400 - dt1.tsec;
            result += (dt2.julnum - 1 - dt1.julnum) * 86400;
        }
        else
        {
            result = dt2.tsec - dt1.tsec;
            result += (dt2.julnum - dt1.julnum) * 86400;
        }
    }
    return result;
} // end operator-


DateTime operator- (const DateTime &dt, int sec)
{
    if (sec < 0)
    {
        return ::operator+(dt, abs(sec));
    }
    else
    {
        if (dt.tsec < (totalSecondsType)sec)
        {
            julType j;
            totalSecondsType ts;
            int needed = (sec - dt.tsec - 1) / (86400) + 1;
            j = dt.julnum - needed;
            ts = dt.tsec + (needed * 86400) - sec;
            return DateTime(j, ts);
        }
        else
        {
            return DateTime(dt.julnum, dt.tsec - sec);
        }
    }
}

DateTime operator- (int sec, const DateTime &dt)
{
    return ::operator-(dt, sec);
}

void operator-= (DateTime &dt, int sec)
{
    if (sec < 0)
    {
        dt += abs(sec);
    }
    else
    {
        if (dt.tsec < (totalSecondsType)sec)
        {
            int needed = ((sec - dt.tsec - 1) / (86400) + 1);
            dt.julnum -= needed;
            dt.tsec = dt.tsec + (needed * 86400) - sec;
        }
        else
        {
            dt.tsec -= sec;
        }
    }
}

int operator<(const DateTime &dt1, const DateTime &dt2)
{
    return int(dt1.julnum < dt2.julnum ||
               (dt1.julnum == dt2.julnum && dt1.tsec < dt2.tsec));
}

int operator<=(const DateTime &dt1, const DateTime &dt2)
{
    return int(dt1 < dt2 || dt1 == dt2);
}


int DateTime::isNull() const
{
    return (*this == NoDateTime);
}

ostream& operator<<(ostream &os, const DateTime &dt)
{
    if (dt.isNull())
    {
        os << "NULL";
    }
    else
    {
        os << (Date&)dt;
        os << (CTime&)dt;
    }
    return os;
}

void DateTime::Print()
{
    char *pt = syCh();
    cout << pt << endl;
    delete pt;
}


void DateTime::getDateTime(yearType &yyyy, monthType &mm, dayType &dd,
                           hourType &hh, minuteType &mn, secondType &ss)
{
    jday2date(dd, mm, yyyy);
    totSec2time(hh, mn, ss);
}

char *DateTime::getChDateTime() const
{
    static char datetime[21];

    sprintf(datetime, "%s", this->getChDate());
    if (strcmp(datetime, "null") != 0)
    {
        sprintf(datetime + 11, " %s", this->getChTime());
    }

    return datetime;
}

char *DateTime::syCh(int form) const
{
    char *pt = new char[22];


    if (*this == NoDateTime)
    {
        sprintf(pt, "null");
    }
    else
    {
        static yearType y = 0;
        static monthType m = 0;
        static dayType d = 0;
        static hourType h = 0;
        static minuteType min = 0;
        static secondType s = 0;
        jday2date(d, m, y);
        totSec2time(h, min, s);
        switch (form)
        {
        case 1:
            sprintf(pt, "%02d/%02d/%d %02d:%02d", d, m, y, h, min);
            break;
        case 2:
            sprintf(pt, "%d%02d%02d%02d%02d", y, m, d, h, min);
            break;
        case 3:
            sprintf(pt, "%d-%02d/%02d %02d:%02d", y, m, d, h, min);
            break;
        case 4:
            sprintf(pt, "%02d/%02d/%d", d, m, y);
            break;
        case 5:
            sprintf(pt, "%04d%02d%02d/%02d%02d", y, m, d, h, min);
            break;
        case 6:
            sprintf(pt, "%02d/%02d %02d:%02d", d, m, h, min);
            break;
        default:
            sprintf(pt, "%d/%d/%d %d:%d:%d", m, d, y, h, min, s);
            break;
        };
    }
    return (pt);
}

void DateTime::internal(char str[13])
{
    static yearType y = 0;
    static monthType m = 0;
    static dayType d = 0;
    static hourType h = 0;
    static minuteType min = 0;
    static secondType s = 0;

    if (isNull())
    {
        sprintf(str, "%s", "null");
    }
    else
    {
        jday2date(d, m, y);
        totSec2time(h, min, s);
        sprintf(str, "%d%02d%02d%02d%02d", y, m, d, h, min);
    }
}

DateTime DateTime::StartOfDay() const
{
    if (*this == NoDateTime)
    {
        return NoDateTime;
    }
    else
    {
        static yearType y = 0;
        static monthType m = 0;
        static dayType d = 0;
        jday2date(d, m, y);
        DateTime dt(y, m, d);
        return dt;
    }
}

DateTime DateTime::StartOfMonth() const
{
    if (*this == NoDateTime)
    {
        return NoDateTime;
    }
    else
    {
        static yearType y = 0;
        static monthType m = 0;
        static dayType d = 0;
        jday2date(d, m, y);
        DateTime dt(y, m, (dayType)1);
        return dt;
    }
}

DateTime DateTime::StartOfYear() const
{
    if (*this == NoDateTime)
    {
        return NoDateTime;
    }
    else
    {
        static yearType y = 0;
        static monthType m = 0;
        static dayType d = 0;
        jday2date(d, m, y);
        DateTime dt(y, (monthType)1, (dayType)1);
        return dt;
    }
}

long int DateTime::difference_minutes(DateTime &dt2)
{
    long int daydiff = date2jday() - dt2.date2jday();
    long int hourdiff = getHour() - dt2.getHour();
    long int mindiff = getMinute() - dt2.getMinute();
    long int diff = daydiff * 24l * 60l + hourdiff * 60l + mindiff;

    return diff;
}

DateTime  &DateTime::add_minutes(long int num_minutes)
{
    totalSecondsType secs = 60 * ((tsec / 60 + num_minutes) % 1440);
    julType days = julnum + (tsec / 60 + num_minutes) / 1440;

    DateTime dt(days, secs);

    *this = dt;

    return *this;
}



DateTime DateTime::EndOfMonth() const
{
    if (*this == NoDateTime)
    {
        return NoDateTime;
    }
    else
    {
        static yearType y = 0;
        static monthType m = 0;
        static dayType d = 0;
        jday2date(d, m, y);
        if ((int)m == 12)
        {
            return EndOfYear();
        }

        DateTime dt(y, m + 1, 1);
        dt -= 60;
        return dt;
    }
}

DateTime DateTime::EndOfDay() const
{
    if (*this == NoDateTime)
    {
        return NoDateTime;
    }
    else
    {
        static yearType y = 0;
        static monthType m = 0;
        static dayType d = 0;
        jday2date(d, m, y);
        DateTime dt(y, m, d, 23, 59);
        return dt;
    }
}

DateTime DateTime::EndOfYear() const
{
    if (*this == NoDateTime)
    {
        return NoDateTime;
    }
    else
    {
        static yearType y = 0;
        static monthType m = 0;
        static dayType d = 0;
        jday2date(d, m, y);
        DateTime dt(y, 12, 31, 23, 59);
        return dt;
    }
}

double DateTime::as_floating_point_year(void)
{
    DateTime start = (*this).StartOfYear();
    long int sec_from_start = (*this) - start;
    long int sec_in_year = 86400 * daysInYear();

    double year = getYear();

    year += ((double)sec_from_start) / ((double)sec_in_year);

    return year;
}

// ******************************************************************
// Class: Datetimes
// take care of converting operations regarding the special DateTime format
// ******************************************************************
const char* Datetimes::Dt2Str(DateTime dt, int type)
{
    // Converts from DateTime to representative string
    // Type represent which format chosen:
    // 1: yyyy/mm/dd hh:mnmn  ex. 1994/05/05 10:55
    // 2: yyyy/mm/dd  ex. 1994/05/05
    // 3: dd/mm/yyyy  ex. 05/05/1994
    // 4: month year  ex. May 1994
    // 5: date month year  ex. 5. May 1994
    // 6: dd/mm/yy hh:mnmn  ex. 13/05/1994 12:05
    // 7: yyyymmdd hhmm   ex. 19951224 1231,  spreadsheet

    if (dt == NoDateTime)
    {
        return("----");
    }

    strcpy(timestr, "");

    switch (type)
    {
    case 1:
        sprintf(timestr, "%4d/%2d/%2d %2d:%2d", dt.getYear(), dt.getMonth(),
                dt.getDay(), dt.getHour(), dt.getMinute());
        break;
    case 2:
        sprintf(timestr, "%4d/%2d/%2d", dt.getYear(), dt.getMonth(), dt.getDay());
        break;
    case 3:
        sprintf(timestr, "%2d/%2d/%4d", dt.getDay(), dt.getMonth(), dt.getYear());
        break;
    case 4:
        sprintf(timestr, "%s %4d", Month[dt.getMonth() - 1], dt.getYear());
        break;
    case 5:
        sprintf(timestr, "%2d. %s %4d", dt.getDay(), Month[dt.getMonth() - 1], dt.getYear());
        break;
    case 6:
        sprintf(timestr, "%2d/%2d/%4d %2d:%2d", dt.getDay(), dt.getMonth(),
                dt.getYear(), dt.getHour(), dt.getMinute());
        break;
    case 7:
        sprintf(timestr, "%04d%02d%02d %02d%02d", dt.getYear(), dt.getMonth(), dt.getDay(),
                dt.getHour(), dt.getMinute());
        break;
    };

    return(timestr);
}


DateTime Datetimes::Str2Dt(char *dtstr)
{
    // Converts a timestring to datetime format
    // input format: 12/11/1994 12:05

    char year[5], month[3], day[3], hour[3], minute[3];

    strncpy(day, dtstr, 2);
    strncpy(month, dtstr + 3, 2);
    strncpy(year, dtstr + 6, 4);
    strncpy(hour, dtstr + 11, 2);
    strncpy(minute, dtstr + 14, 2);

    day[2] = '\0';
    month[2] = '\0';
    year[4] = '\0';
    hour[2] = '\0';
    minute[2] = '\0';

    if (!strlen(dtstr))
    {
        return NoDateTime;
    }

    return(*new DateTime(atoi(year), atoi(month), atoi(day),
                         atoi(hour), atoi(minute), 0));

}


char* DateTime::getRfc822Str(int UT_offset) const
{
    const char *month_names_eng[12] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
                                };
    char* strPtr = new char[32];

    sprintf(strPtr, "%02d %3s %04d %02d:%02d:%02d %c%02d00",
            getDay(), month_names_eng[getMonth() - 1], getYear(),
            getHour(), getMinute(), getSecond(),
            UT_offset >= 0 ? '+' : '-', abs(UT_offset));

    return strPtr;
}

int overlap(const DateTime &start1, const DateTime &end1,
            const DateTime &start2, const DateTime &end2)
{
    if (start1 != NoDateTime && end2 != NoDateTime && start1 > end2)
    {
        return 0;
    }

    if (end1 != NoDateTime && start2 != NoDateTime && end1 < start2)
    {
        return 0;
    }

    return 1;
}


#ifdef MAIN
main()
{
    DateTime d(1953, 1, 1), dd;
    cout << d << endl;
    cout << dd << endl;
    cout << d.getRfc822Str(-3) << endl;

    DateTime a(2004, 4, 5, 4, 23);
    int pos = 46000;
    int neg = -46000;
    cout << "a + 46000 = " << a + 46000 << endl;
    cout << "46000 + a = " << 46000 + a << endl;
    cout << "a - 46000 = " << a - 46000 << endl;
    cout << "46000 - a = " << 46000 - a << endl;
    cout << "a - -46000 = " << a - -46000 << endl;
    cout << "-46000 - a = " << -46000 - a << endl;
    cout << "a + -46000 = " << a + -46000 << endl;
    cout << "-46000 + a = " << -46000 + a << endl;

}
#endif
