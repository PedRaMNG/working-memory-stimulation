function sendTelegramMessage(message)
   % Your personal Telegram credentials (ask ChatGPT to give you instructions on setting it up)
   botToken = ''; %should be changed
   chatID = ''; %should be changed
   
   % Create the URL
   url = sprintf('https://api.telegram.org/bot%s/sendMessage', botToken);
   
   % Prepare the data
   data = struct('chat_id', chatID, 'text', message);
   
   % Convert data to JSON
   options = weboptions('MediaType', 'application/json');
   
   try
       % Send the HTTP request
       response = webwrite(url, data, options);
       fprintf('Message sent successfully!\n');
   catch exception
       fprintf('Failed to send message: %s\n', exception.message);
   end
end

% Example usage:
% sendTelegramMessage('Your MATLAB simulation has completed!');