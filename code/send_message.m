function send_message(id,subject,message,attachment)
% SEND_MESSAGE send email or text message
%    SEND_TEXT_MESSAGE(ID,SUBJECT,MESSAGE,ATTACHMENT/CARRIER) 
%       ID -  Email address of 10-digit cell phone number.
%       ATTACHMENT/CARRIER - Either filename(s) of attachment(s) or cell
%           phone service provider:'Alltel', 'AT&T', 'Boost',
%           'Nextel', 'Sprint', 'T-Mobile', 'Verizon', or 'Virgin'.
%
%   See also SENDMAIL.

mail = 'TPFLab@gmail.com';
password = 'TenToTheTen';

if nargin < 2
    error('At least two inputs are needed');
elseif nargin == 2
    message = [];
    attachment = [];
elseif nargin == 3
    attachment = [];
end

if ~isempty(attachment)
    sms = true;
    switch strrep(strrep(lower(attachment),'-',''),'&','')
        case 'alltel',    id = strcat(strrep(id, '-', ''),'@message.alltel.com');
        case 'att',       id = strcat(strrep(id, '-', ''),'@txt.att.net');
        case 'boost',     id = strcat(strrep(id, '-', ''),'@myboostmobile.com');
        case 'nextel',    id = strcat(strrep(id, '-', ''),'@messaging.nextel.com');
        case 'sprint',    id = strcat(strrep(id, '-', ''),'@messaging.sprintpcs.com');
        case 'tmobile',   id = strcat(strrep(id, '-', ''),'@tmomail.net');
        case 'verizon',   id = strcat(strrep(id, '-', ''),'@vtext.com');
        case 'virgin',    id = strcat(strrep(id, '-', ''),'@vmobl.com');
        otherwise,        sms = false;
    end
    if sms, attachment = []; end
end

% Set up Gmail SMTP service.
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);

% Gmail server.
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail(id,subject,message,attachment)
end
